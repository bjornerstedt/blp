classdef RCDemand < NLDemand
    % Random Coefficient demand class
    %   Simulation of shares and estimation of parameters beta and sigma.
   
    properties (SetAccess = protected, Hidden = true )
        draws % Random draws (nind, nmkt, k) matrix k = # random params. 
        W
        x2 % nonlinear parameter vector
        xi % Unobservable utility, residual
        vars2 % Nonlinear variable names in output (with rc_ prefix)
        nonlinparams = [] % Arrays for RC names and type of RC
        nonlintype = [] % Array for type of RC
        oldsigma = 0 % Used in comparisons between minimization steps
        deltaJac % Used to guess new edelta
        isOptimalIV = false
        edelta % Saved between invocations of objective
        v
        iweight % quadrature weights     
        ZWZ
        inv_x1ZWZx1
        estimationMatrix
    end
    
    methods
%% GENERAL DEMAND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function [s, si] = shares(obj, p)
            [s, si] = obj.period{obj.sim.market}.shares( p );
        end
       
        function [sh, indsh] = sharesAll(obj, p)
            sel = logical(obj.dummarket);
            sh = zeros(size(p));
            indsh = zeros(length(p), obj.settings.nind);
            for t = 1:max(obj.marketid)
                [s, si] = obj.period{t}.shares( p(sel(:,t)) );
                sh(sel(:,t)) = s;
                indsh(sel(:,t)) = si;
            end
        end
       
		% Calculate quantities and market shares
        function tab = quantity(obj, P)
            tableCols = {'Product', 'Price', 'Quantity', 'MarketSh', 'Share'};
            s = obj.shares(P);
            if obj.settings.ces
                q = s .* obj.ms ./ P;
            else
                q = s .* obj.ms;
            end
            m = q / sum(q); 
            tab = table(obj.panelid, P, q, m, s);
            tab.Properties.VariableNames = tableCols;
        end
                
        function sh = shareJacobian(obj, P)
            % shareJacobian(p, [market_number])
            sh = obj.period{obj.sim.market}.shareJacobian(P);          
        end
        
        function [elas] = groupElasticities(obj, P,  group)
            group = obj.data{:, group};
            if iscategorical(group)
                [names,~,group] = unique(group);
                names = matlab.lang.makeValidName(cellstr(char(names)));
            else
                names = [];
            end
            A = dummyvar(group)';
            s = obj.shares(P);
            D = obj.shareJacobian(P)';
            elas = A*diag(P)* D * A' * diag( 1 ./(A*s) ) ;
            elas = array2table(elas);
            if ~isempty(names)
                elas.Properties.RowNames = names;
                elas.Properties.VariableNames = names;
            end
        end
        
        function s = sumstats(obj, j, E)
            e = j .* E;
            e(e==0)=[];
            s = [mean(e) std(e) min(e) max(e)];
        end
        
        function [elas, varargout] = elasticities(obj, P, varargin)
            s = obj.shares(P);
            n = length(s);
            D = obj.shareJacobian(P)';
            E = diag(P) * D * diag( 1 ./ s );
            elas = obj.sumstats(eye(n), E);
            if nargin == 0
                elas = [elas; obj.sumstats(1 - eye(n), E)];
                rowtit = {'e_ii', 'e_ij'};
            else
                G = dummyvar( obj.data{:, varargin{1}} );
                GG = G*G';
                elas = [elas; ...
                    obj.sumstats(GG - eye(n), E); ...
                    obj.sumstats(1 - GG, E)];
                rowtit = {'e_ii', 'e_ij', 'e_ik'};
            end
            elas = array2table(elas);
            elas.Properties.VariableNames = {'Mean', 'Std', 'Min', 'Max'};
            elas.Properties.RowNames = rowtit;
            if nargout > 0
                varargout{1} = E;
            end
        end
        
        function der = deltaJacobian(obj, sigma, edelta)
            der = zeros(size(obj.x2));
            for t = 1:max(obj.marketid)
                index = obj.marketid == t;
                der(index, :) = obj.period{t}.deltaJacobian(sigma, edelta);
            end
        end
        
%% Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function sigma = minimize(obj, varargin)
        %        extraoptions = {'OutputFcn', @iterOutput}; % Show output func
% 'DerivativeCheck', 'on', ...
%             obj.config.historyx = [];
%             obj.config.historyfval = [];
            if nargin > 1
                extraoptions = varargin{1};
            else
                extraoptions = {};
            end
            if isempty(obj.sigma)
                error('sigma is not set, probably because init has not been invoked');
            end
            obj.oldsigma = zeros(size(obj.sigma));
                options = optimoptions(@fminunc, ...
                    'MaxIter',obj.settings.maxiter, ... 
                    'TolX', obj.config.tolerance, 'TolFun', obj.config.tolerance, ...
                    'Algorithm','trust-region' ,...
                    'GradObj','on', extraoptions{:});
            fval = 10^6;
            i = 0;
            finished = false;
            func = @(x)obj.objective(x);
            while ~finished && i <= obj.config.restartMaxIterations 
                [sigma,fval,exitflag] = fminunc(func, obj.sigma, options);
                if obj.settings.optimalIV && fval > obj.config.restartFval || exitflag <= 0
                    obj.get_sigma();
                    if exitflag <= 0
                        disp('Restarting as minimization did not converge.');
                    else
                        disp(['Restarting as fval=', num2str(fval), ' is too large']);
                    end
                else
                    finished = true;      
                end
                i = i + 1;
            end
            obj.results.other.fval = fval;
            obj.results.other.restarts = i - 1;
            obj.results.other.minimum = finished;
        end
        
        function [f, g] = objective(obj, sigma)
            %This function defines the objective over which to minimize
            obj.edelta =  obj.findDelta(sigma);
            del = log(obj.edelta);
            
            if max(isnan(del)) == 1
                f = 1e+10;
                if nargout > 1
                    g = 1e+10;
                end
            else
                [~, xi] = obj.lpart(del);
                if ~isempty(obj.Z)
                    f = xi' * obj.ZWZ * xi;
                    if nargout > 1
                        obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
                        g = 2* obj.deltaJac'* obj.ZWZ * xi;
                    end
                else
%                     xiX =  xi' * obj.X;
%                     f = xiX * xiX';
                    f = xi' * xi;
                    if nargout > 1
                         obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
%                         g = 2*obj.deltaJac'* obj.X * xiX';
                        g = 2*obj.deltaJac' * xi;
                    end
                end
            end
        end
        
        function R = estimate(obj, varargin)
            obj.init(varargin{:});
            obj.edelta = exp(obj.share.ls);
            est = obj.estimation_step(false, varargin{:});
            if obj.settings.optimalIV
                R = obj.estimation_step(true, varargin{:});
                obj.results.estimate1 = est;   
            else
                R = est;
            end
       end
        
        function R = estimation_step(obj, optIV, varargin)
            obj.initEstimation(optIV);
            obj.sigma = obj.minimize([{'Display','iter-detailed'}, varargin]);
            delta = log(obj.edelta);
            if isnan(delta)
                error('delta is NaN');
            end
            obj.sigma = abs(obj.sigma); % Works for symmetric or positive dist
            [obj.beta, xi] = obj.lpart(delta);
            
            if strcmpi(obj.settings.paneltype, 'fe')
                % xi is the residual of a regression of demeaned vars on
                % delta, which is not demeaned
                davt = accumarray(obj.panelid, delta, [], @mean);
                davt = davt(obj.panelid, :);
                obj.xi = xi - davt;
            else
                obj.xi = xi;
            end            
            
            obj.alpha = - obj.beta(strcmp(obj.getPriceName(), obj.vars));
            coef = [obj.beta; obj.sigma];  % vector of all parameters
            varsel = [1:length(obj.vars ), ...
                (length(coef)-length(obj.vars2)+1):length(coef)];
            varcovar = obj.computeVariance();
            varcovar = varcovar(varsel, varsel);
            variables = [obj.vars obj.vars2];
            obj.results.estimate = obj.createResults(coef(varsel), varcovar, variables);
            if strcmpi(obj.settings.paneltype, 'lsdv')
                obj.results.betadummies = ...
                    coef((length(obj.vars )+1):(length(coef)-length(obj.vars2)));
            end
            varcovar = array2table(varcovar );
            varcovar.Properties.VariableNames = variables;
            varcovar.Properties.RowNames = variables;
            obj.results.params.varcovar = varcovar;
            obj.results.params.sigma = obj.sigma;
            R = obj.results.estimate;
        end
        
        function varcovar = computeVariance(obj)
            derdel = obj.deltaJacobian(obj.sigma, obj.edelta);
            if strcmpi(obj.settings.paneltype, 'fe')
                dgf = (size(obj.X,1) - size(obj.X,2)) - max(obj.panelid);
            else
                dgf = (size(obj.X,1) - size(obj.X,2));
            end
            obj.results.params.dgf = dgf;
            if isempty(obj.Z)
                % Homoscedastic errors
                dertheta = [-obj.X derdel];
                invXX = inv(dertheta' * dertheta);
                if obj.settings.robust
                    xiX =bsxfun(@times, obj.xi, dertheta);
                    varcovar = invXX*(xiX'*xiX)*invXX;
                else
                    varcovar = ((obj.xi'*obj.xi) ./ dgf) * invXX;
                end
            else
                if obj.config.test
                    pHat = obj.Z*((obj.Z'*obj.Z)\obj.Z'*obj.X(:, 1));
                    deltaHat = [pHat, obj.X(:,2:end)] * obj.beta;
                end
                dertheta = [-obj.X derdel]'* obj.Z;
% Robust estimation if test is set:
%             if isfield(obj.config, 'test') && obj.config.test
%                 xiZ = bsxfun(@times, obj.xi, obj.Z );
%                 obj.W = inv(xiZ'*xiZ);
%             end
                obj.results.other.cond = cond(dertheta * obj.W * dertheta');
                if obj.settings.robust
                    varcovar = inv(dertheta * obj.W * dertheta');
                else
                    varcovar = (obj.xi' * obj.xi/dgf) .* inv(dertheta * obj.W * dertheta');
                end
            end
        end
                      
%% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function init(obj, varargin)
            % init(data, selection) - both arguments optional    
            selection = init@NLDemand(obj, varargin{:});
            obj.nonlinparams = [];
            obj.nonlintype = [];
            nonlin = {obj.var.nonlinear, obj.var.nonlinearlogs, ...
                obj.var.nonlineartriangular};
%             for i = 1:length(nonlin)
            for i = 1:1 % HACK !!!!!!!!!!!!!!!!
                nvars = strtrim(nonlin{i});
                if ~isempty(nvars)
                    params = strsplit(nvars);
                    obj.nonlintype = ones(length(params),1)*(i - 1);
                    obj.nonlinparams = [obj.nonlinparams, params];
                end
            end
            if isempty(obj.nonlinparams)
                error('Some variable has to be specified as nonlinear');
            end
            if isempty(obj.sigma)
                obj.get_sigma();
            end
            if size(obj.sigma, 2) > 1
                obj.sigma = obj.sigma';
            end
            if length(obj.nonlinparams) ~= length(obj.sigma)
                error('sigma and nonlinear have different dimensions');
            end
            obj.results.sigma0 = obj.sigma;
            obj.vars2 = ...
                cellfun(@(x) {sprintf('rc_%s', x)}, obj.nonlinparams );
            obj.x2 = obj.data{:, obj.nonlinparams };
            nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
            if any(nonlinprice) 
                % For CES, logprice rather than price.
                obj.x2(:, nonlinprice) = obj.X(:, 1);
            end
            obj.randdraws(selection);
            if ~isempty(selection)
                obj.x2 = obj.x2(selection, :); 
            end
            obj.initPeriods();
        end
        
        function get_sigma(obj)
            if isempty(obj.config.randstream)
                obj.sigma = randn(length(obj.nonlinparams),1);
            else
                obj.sigma = obj.config.randstream.randn(length(obj.nonlinparams),1);
            end
        end
        
        function initPeriods(obj)
            md = RCDemandMarket(obj);
            obj.period = cell(max(obj.marketid), 1);
            for t = 1:max(obj.marketid)
                newmarket = copy(md);
                newmarket.selection = logical(obj.dummarket(:, t));
                newmarket.x2 = obj.x2(newmarket.selection, :);
                if ~isempty(obj.share)
                    newmarket.s = obj.share.s(newmarket.selection, 1);
                end
                if ~isempty(obj.d)
                    newmarket.d = obj.d(newmarket.selection, 1);
                end
                newmarket.p = obj.p(newmarket.selection, 1);
                newmarket.v = [];
                for k = 1:size(obj.x2, 2)
                    temp = obj.v(newmarket.selection,:,k);
                    newmarket.v = temp(1,:); 
                end
                newmarket.init();
                obj.period{t} = newmarket;
            end
        end
            
        function randdraws(obj, selection)
            K = size(obj.x2, 2); % Number of variables
            obj.v = [];
            if strcmpi(obj.settings.drawmethod, 'quadrature')
                [X, obj.iweight] = nwspgr('KPN', K, obj.settings.quaddraws);
                obj.settings.nind = length(obj.iweight);
                for k = 1:K
                    obj.v(:,:,k) = repmat(X(:,k)', size(obj.x2, 1), 1);
                end
            else
                obj.iweight = ones(obj.settings.nind, 1) / obj.settings.nind ;
                if isempty(obj.draws) % Can be set manually
                    if obj.settings.marketdraws
                        obj.draws = Sampling.draw(obj.settings.drawmethod, ...
                            K, max(obj.marketid) * obj.settings.nind, ...
                            obj.config.randstream);
                    else
                        obj.draws = Sampling.draw(obj.settings.drawmethod, ...
                            K, obj.settings.nind, obj.config.randstream);
                    end
                end
                for k = 1:K
                    % For each k, create draws for each market - wk (nind x nmkt)
                    if obj.settings.marketdraws
                        wk = reshape( obj.draws(:, k), obj.settings.nind, ...
                            max(obj.marketid))';
                    else
                        wk = repmat( obj.draws(:, k)', max(obj.marketid), 1);
                    end
                    % Duplicate draws for all products in each market
                    % Each k can have a different type of distribution
                    if obj.nonlintype(k) == 0
                        obj.v(:,:,k) = wk(obj.marketid, :);
                    elseif obj.nonlintype(k) == 1 % log-normal
                        obj.v(:,:,k) = exp(wk(obj.marketid, :)); 
                    elseif obj.nonlintype(k) == 2 % triangular
                        obj.v(:,:,k) = Sampling.triangular(wk(obj.marketid, :)); 
                    end              
                end
            end
            if ~isempty(selection)
                vnew = [];
                % This could be done above in creation, but this is simpler
                for k = 1:K
                    vnew(:,:,k) = obj.v(selection, :, k);
                end
                obj.v = vnew;
            end
        end
        
        function initSimulation(obj, market)
            % General init, move to estimate? Should NestedLogit have
            % similar code?
            obj.sim.market = market; % HACK to get market based routines to work
            if isempty(obj.d)
                % Create starting values for findDelta
                obj.edelta = obj.findDelta(obj.sigma);
                obj.d = log(obj.edelta) + obj.alpha*obj.X(:, 1);
                for t = 1:max(obj.marketid)
                    obj.period{t}.d = obj.d(obj.dummarket(:, t));
                end
            end
        end
                                  
        function demand = pack(obj, filename)
            demand = RCDemand;
            demand.settings = obj.settings;
            demand.var = obj.var;
            demand.config = obj.config;
            demand.results = obj.results;
            demand.beta = obj.beta;
            demand.sigma = obj.sigma;
        end
        
        function obj = RCDemand(varargin)
            obj = obj@NLDemand(varargin{:});
            obj.var.setParameters({'nonlinear','nonlinearlogs','nonlineartriangular'});
            obj.settings.setParameters({'optimalIV','drawmethod',... 
                'marketdraws','nind','quaddraws','maxiter','fptolerance1','fptolerance2'});
            obj.config = SettingsClass({'tolerance','fpmaxit', ...
                'restartMaxIterations','restartFval', 'test', ...
                'randstream','hessian','guessdelta','quietly'});
            
            obj.settings.paneltype = 'lsdv';

            obj.settings.optimalIV = false;
            obj.settings.drawmethod = 'hypercube';
            obj.settings.marketdraws = false; % Different draws for each market
            obj.settings.nind = 100; % number of simulated "indviduals" per market 
            obj.settings.quaddraws = 10;
            obj.settings.maxiter = 100;
            obj.settings.fptolerance1 = 1e-14; % lower tol for first FP its
            obj.settings.fptolerance2 = 1e-14; % high tol for last iterations
            obj.settings.ces = false;
            
            obj.config.tolerance = 1e-9;
            obj.config.fpmaxit = 1000; % maximum iterations in contraction mapping
            obj.config.restartMaxIterations = 1; % Max # of restarts if fval is high
            obj.config.restartFval =  10^3; % Restart optimal IV estimation if fval >.
            obj.config.randstream = []; % random stream for estimation in parallel 
            obj.config.guessdelta = true;
            obj.config.hessian = false;
            obj.config.quietly = true;
            obj.results.estimateDescription = 'Random Coefficient Logit Demand'; 
        end  
        
    end

%% Basic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    methods % (Access = private)
        function newedelta = findDelta(obj, sigma)
            if isempty(obj.edelta)
                % Create starting values for findDelta
                if isempty(obj.share.s)
                    obj.share.s = sum(bsxfun(@times, obj.dummarket, ...
                        1 ./ sum(obj.dummarket) ), 2) .* 0.5;
                    obj.edelta = obj.share.s .* 2;
                else
                    obj.edelta = obj.share.s ./ obj.share.s0;
                end
            end
            if max(abs(sigma - obj.oldsigma)) < 0.01
                tolerance = obj.settings.fptolerance2;
                closeFlag = 0;
            else
                tolerance = obj.settings.fptolerance1;
                closeFlag = 1;
            end
            
            if obj.config.guessdelta && ~isempty(obj.deltaJac)
                newdelta = log(obj.edelta) + ...
                    obj.deltaJac*(sigma - obj.oldsigma); 
                edelta = exp(newdelta);
            else
                edelta =  obj.edelta;
            end

            sel = logical(obj.dummarket);
            newedelta = zeros(size(edelta));
            for t = 1:max(obj.marketid)
                newedelta(sel(:,t)) = obj.period{t}.findDelta(sigma, ...
                    edelta(sel(:,t)), tolerance);
            end
            % Update oldsigma and edelta only in first stage and if
            % successful
            if closeFlag == 1 && max(isnan(newedelta)) < 1;
                obj.oldsigma = sigma;
            end
        end       
        
        function [ bet, xi ] = lpart(obj, del )
%             if strcmpi(obj.settings.paneltype, 'fe')
%                 davt = accumarray(obj.panelid, del, [], @mean);
%                 davt = davt(obj.panelid, :);
%                 del = del - davt;
%             end
            if ~isempty(obj.Z)
                bet = obj.inv_x1ZWZx1 * (obj.X' * obj.ZWZ * del);
            else
                bet = obj.inv_x1ZWZx1 * (obj.X' * del);
            end
%             bet = obj.estimationMatrix * del;

            xi = del - obj.X * bet;
        end
        
        function initEstimation(obj, optIV)
            if ~optIV
                if isempty(obj.W) && ~isempty(obj.Z) % W can be specified manually
                    obj.W = inv(obj.Z' * obj.Z);
                end
                obj.isOptimalIV = false;
            else
                if obj.isOptimalIV
                    return % Already calculated optimal IV
                end
                obj.isOptimalIV = true;
                %                 pHat = obj.Z*((obj.Z'*obj.Z)\obj.Z'*obj.X(:, 1));
                %                 deltaHat = [pHat, obj.X(:,2:end)] * obj.beta;
                % The following code is an attempt to get optimal IV with
                % FE to give the same result as LSDV
                pHat = obj.Zorig*((obj.Zorig'*obj.Zorig)\obj.Zorig'*obj.Xorig(:, 1));
                deltaHat = [pHat, obj.Xorig(:, 2:end)] * obj.beta;
                optInstr = obj.deltaJacobian(obj.sigma, exp(deltaHat));
                if strcmpi(obj.settings.paneltype, 'fe')
                    da = [];
                    for i = 1:size(optInstr ,2)
                        da = [da, accumarray(obj.panelid, optInstr(:, i),...
                            [],@mean)];
                    end
                    davt = da(obj.panelid, :);
                    optInstr  = (optInstr  - davt);
                end                
%                 phatMean = accumarray(obj.panelid, pHat,[],@mean);
%                 pHat = pHat  - phatMean(obj.panelid, :);
                obj.Z  = [obj.X(:,2:end), pHat, optInstr];
                xiZ = bsxfun(@times, obj.xi, obj.Z );
                obj.W = inv(xiZ'*xiZ);
            end
            
%             if ~isempty(obj.Z)
%                 obj.ZWZ = obj.Z * obj.W * obj.Z';
%                 obj.estimationMatrix = inv(obj.X'*obj.ZWZ*obj.X) * obj.X' * obj.ZWZ;
%             else
%                 obj.estimationMatrix = inv(obj.X'*obj.X) * obj.X';
%             end
            if ~isempty(obj.Z)
                obj.ZWZ = obj.Z * obj.W * obj.Z';
                obj.inv_x1ZWZx1 = inv(obj.X'*obj.ZWZ*obj.X);
            else
                obj.inv_x1ZWZx1 = inv(obj.X'*obj.X);
            end
        end
        
    end
end

