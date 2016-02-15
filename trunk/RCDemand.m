classdef RCDemand < NLDemand
    % Random Coefficient demand class
    %   Simulation of shares and estimation of parameters beta and sigma.
   
%     properties (SetAccess = protected, Hidden = true )
    properties (SetAccess = protected )
        draws % Random draws (nind, nmkt, k) matrix k = # random params. 
        W
        x2 % nonlinear parameter vector
        xi % Unobservable utility, residual
        vars2 % Nonlinear variable names in output (with rc_ prefix)
        nonlinparams = [] % Arrays for RC names and type of RC
        
        edelta % Saved between invocations of objective
        oldsigma = 0 % Used in comparisons between minimization steps
        deltaJac % Used to guess new edelta
        isOptimalIV = false
        ZWZ
        inv_x1ZWZx1
%         estimationMatrix
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
            if isempty(P) 
                P = obj.p(obj.dummarket(:, obj.sim.market));
            end
            sh = obj.period{obj.sim.market}.shareJacobian(P);          
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
                    obj.getSigma();
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
                disp('Negative deltas found at observations')
                err = find(isnan(del));
                if length(err) < 11
                    disp(err)
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
            est = obj.estimationStep(false, varargin{:});
            if obj.settings.optimalIV
                R = obj.estimationStep(true, varargin{:});
                obj.results.estimate1 = est;   
            else
                R = est;
            end
       end
        
        function R = estimationStep(obj, optIV, varargin)
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
            % randdraws can only be run once, in order to not get new draws
            % in simulation after estimation
            if isempty(obj.draws)
                obj.randdraws();
            end
            obj.x2 = obj.data{:, obj.nonlinparams };
            % log price in CES in x2:
            nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
            if any(nonlinprice) && obj.settings.ces
                obj.x2(:, nonlinprice) = log(obj.data{:, obj.var.price});
                nonlinparamsCES = obj.nonlinparams;
                nonlinparamsCES{nonlinprice} = 'lP';
                obj.vars2 = ...
                    cellfun(@(x) {sprintf('rc_%s', x)}, nonlinparamsCES );
            else
                obj.vars2 = ...
                    cellfun(@(x) {sprintf('rc_%s', x)}, obj.nonlinparams );
            end
            if ~isempty(selection)
                obj.x2 = obj.x2(selection, :); 
            end
            obj.initPeriods();
        end
        
        function randdraws(obj)
            if obj.settings.marketdraws
                markets = max(obj.marketid);
            else
                markets = 1;
            end
            obj.draws = Draws('DrawMethod', obj.settings.drawmethod,...
                'Markets', markets, 'Individuals', obj.settings.nind,...
                'Accuracy', obj.settings.quaddraws, ...
                'RandStream', obj.config.randstream);
            if isempty(obj.var.nonlinear)
                error('Some variable has to be specified as nonlinear');
            end
            % The parse method could be included in create, but this will
            % affect draw orders and thus all tests. A little work ...
            obj.nonlinparams = obj.draws.parse( obj.var.nonlinear);
            obj.draws.create( );
            obj.getSigma();
        end
        
        function getSigma(obj)
            % Only run once, but why? Remove in cleaning obj.sigm
            if isempty(obj.sigma)
                if ~isempty(obj.settings.sigma0)
                    obj.sigma = obj.settings.sigma0;
                elseif isempty(obj.config.randstream)
                    obj.sigma = randn(length(obj.nonlinparams), 1);
                else
                    obj.sigma = obj.config.randstream.randn(length(obj.nonlinparams),1);
                end
            end
            obj.results.sigma0 = obj.sigma;
            if size(obj.sigma, 2) > 1
                obj.sigma = obj.sigma';
            end
            if length(obj.nonlinparams) ~= length(obj.sigma)
                error('sigma and nonlinear have different dimensions');
            end
        end
        
        function initPeriods(obj)
            md = RCDemandMarket(obj);
            md.iweight = obj.draws.weights;
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
                newmarket.v = [];
                if obj.settings.marketdraws
                    newmarket.v = obj.draws.draws(:,:,t)'; 
                else
                    newmarket.v = obj.draws.draws'; 
                end
                newmarket.init();
                obj.period{t} = newmarket;
            end
        end
              
        function initSimulation(obj, market)
            % General init, move to estimate? Should NestedLogit have
            % similar code?
            obj.sim.market = market; % HACK to get market based routines to work
            if isempty(obj.d)
                % Create starting values for findDelta
                obj.edelta = obj.findDelta(obj.sigma);
                obj.d = log(obj.edelta) + obj.alpha * obj.Xorig(:, 1);
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
                'marketdraws','nind','quaddraws','maxiter', 'sigma0'});
            obj.config = SettingsClass({'tolerance','fptolerance1','fptolerance2', ...
                'restartMaxIterations','restartFval', 'test', 'fpmaxit',...
                'randstream','hessian','guessdelta','quietly'});
            
            obj.settings.paneltype = 'lsdv';

            obj.settings.optimalIV = false;
            obj.settings.drawmethod = 'hypercube';
            obj.settings.marketdraws = false; % Different draws for each market
            obj.settings.nind = 100; % number of simulated "indviduals" per market 
            obj.settings.quaddraws = 10;
            obj.settings.maxiter = 100;
            obj.settings.ces = false;
            obj.settings.sigma0 = [];
            
            obj.config.tolerance = 1e-9;
            obj.config.fptolerance1 = 1e-14; % lower tol for first FP its
            obj.config.fptolerance2 = 1e-14; % high tol for last iterations
            obj.config.fpmaxit = 1000; % maximum iterations in contraction mapping
            obj.config.restartMaxIterations = 1; % Max # of restarts if fval is high
            obj.config.restartFval =  10^3; % Restart optimal IV estimation if fval >.
            obj.config.randstream = []; % random stream for estimation in parallel 
            obj.config.guessdelta = true;
            obj.config.hessian = false;
            obj.config.quietly = true;
            obj.results.estimateDescription = 'Random Coefficient Logit Demand'; 
            obj.results.sigma0 = []; 
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
                tolerance = obj.config.fptolerance2;
                closeFlag = 0;
            else
                tolerance = obj.config.fptolerance1;
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
                % Note that Xorig has log(p) for CES
                pHat = obj.Zorig*((obj.Zorig'*obj.Zorig)\obj.Zorig' * obj.Xorig(:, 1));
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

