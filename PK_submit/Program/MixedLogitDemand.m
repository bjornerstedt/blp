classdef MixedLogitDemand < DemandEstimate
    % Random Coefficient demand class
    %   Simulation of shares and estimation of parameters beta and rc_sigma.
    %   $Id: MixedLogitDemand.m 121 2015-06-08 08:37:01Z d3687-mb $
   
    properties
        alpha     % beta defined in Estimate    
        rc_sigma  % Starting and end value
        draws % Random draws (nind.nmkt, k) matrix k = # random params. 
    end
    properties % (SetAccess = private )
        xi % Residual utility
        d % Utility without the price effect (Make protected?)
        x0 % x1 without prices
        x2 % nonlinear parameter vector

        W
        inv_x1ZWZx1
        ZWZ
        optimalInstrument = []
        v
        vx
        quadw % quadrature weights
        
        nonlinparams = [] % Arrays for RC names and type of RC
        vars2 % Nonlinear variable names in output (with rc_ prefix)
        nonlintype = [] % Array for type of RC
        
        edelta 
        expmu
        oldsigma = 0 % Used in comparisons between minimization steps
        deltaJac % Used to guess new edelta
    end
    
    methods
%% GENERAL DEMAND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function [s, varargout] = shares(obj, p, varargin)
            if obj.settings.ces
                p = log(p);
            end            
            nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
            if any(nonlinprice)
                x2pred = obj.x2; 
                x2pred(:, nonlinprice) = p;
                for k = 1:size( x2pred, 2)
                    obj.vx(:,:,k) = bsxfun(@times, x2pred(:,k), obj.v(:,:,k));
                end
            end
            obj.nlpart(obj.rc_sigma); 
            if nargin == 2 
                edeltapred = exp(obj.d + obj.alpha*p); % obj.d set in init()
            else
                % Used when simulating without estimate, predicted value as
                % in simulating a market.
                edeltapred = exp([p obj.X(:, 2:end)]*obj.beta + varargin{1}); % predicted shares
            end
            [s, indshare] = obj.share(edeltapred);
            if nargout >1
                varargout{1} = indshare;
            end
        end
        
        function S = actualDemand(obj)
            % Returns actual shares, used in Market.findcosts()
            S = obj.s;	
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
                
        function jac = shareJacobian(obj, P)
            if isempty(P) 
                P = obj.p;
            end
            if obj.config.analyticJacobian
                jac = obj.shareJacobianAnal( P);
            else
                func = @(p)(obj.shares(p));
                jac = jacobian(func, P);
            end
            if obj.settings.ces
                S = obj.shares(P);
                jac = (jac - diag(S) ) * diag(1./P) ;
            end            
        end
        
        function sh = shareJacobianAnal(obj, P)
            [~, si] = obj.shares(P);
            if ~obj.config.quadrature
                obj.quadw = ones(size(si,2), 1) / obj.settings.nind ;
            end            
            obj.nonlinparams = strsplit(strtrim(obj.var.nonlinear));
            nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
            if any(nonlinprice) 
                theta_p = obj.rc_sigma(nonlinprice);
                dUi = (obj.alpha + theta_p * obj.v(:,:,nonlinprice));
                ms = dUi .* si;
                sh = zeros(size(si,1));
                for i = 1:size(si,2)
                    aS = ms(:, i);
                    sh = sh + obj.quadw(i)*( diag(aS) - si(:, i) * aS' );
                end
            else
                sh = zeros(size(si,1));
                for i = 1:size(si,2)
                    dUi = obj.alpha*obj.quadw(i);
                    shi = si(:, i);
                    sh = sh + dUi *( diag(shi) - shi*shi' );
                end
            end
        end
        
        function [elas] = groupElasticities(obj, P,  group)
            group = obj.T{:, group};
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
                G = dummyvar( obj.T{:, varargin{1}} );
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
        
        function d = deltaJacobian(obj, rc_sigma, delta0)
            [~, sh] = obj.share(delta0);
            der = zeros(size(obj.x2));
            for t = 1:max(obj.marketid)
                index = obj.marketid == t;
                si = sh(index, :);
                if obj.config.quadrature
                    wsi =  bsxfun(@times, si , obj.quadw'); 
                else
                    wsi = si ./ obj.settings.nind;
                end
                dssigma = zeros(size(si, 1), length(rc_sigma));
                for k = 1:length(rc_sigma)
                    svx = wsi .* obj.vx(index,:,k); 
                    dssigma(:, k) = sum(svx, 2) - si *sum(svx)';
                end
                dsdelta = diag(sum(wsi, 2)) - wsi*si';
                % Note that dividing by obj.settings.nind cancels out.
                dj = -dsdelta\dssigma;
                der(index, :) = dj;
            end
            d = der;
        end
        
%% Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function fval = minimize(obj, varargin)
            if nargin > 1
                extraoptions = varargin{1};
            else
                extraoptions = {};
            end
            if isempty(obj.rc_sigma)
                error('rc_sigma is not set, probably because init has not been invoked');
            end
            obj.oldsigma = zeros(size(obj.rc_sigma));
            fval = 10^6;
            i = 0;
            finished = false;
            while ~finished && i <= obj.config.restartMaxIterations 
                options = optimoptions(@fminunc, ...
                    'MaxIter',obj.settings.maxiter, ... 
                    'TolX', obj.config.tolerance, 'TolFun', obj.config.tolerance, ...
                    'Algorithm','trust-region' ,...
                    'GradObj','on', extraoptions{:});
                [rc_sigma,fval,exitflag] = fminunc( ...
                    @objective, obj.rc_sigma, options);
                if obj.settings.optimalIV && fval > obj.config.restartFval || exitflag <= 0
                    if isempty(obj.config.randstream)
                        obj.rc_sigma = randn(size(obj.rc_sigma));
                    else
                        obj.rc_sigma = obj.config.randstream.randn(size(obj.rc_sigma));                    
                    end
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

            obj.rc_sigma = rc_sigma;
            
            %%%%% Sub functions
            function [f, g] = objective(rc_sigma)
                %This function defines the objective over which to minimize
                if max(abs(rc_sigma - obj.oldsigma)) < 0.01
                    tolerance = obj.settings.fptolerance2;
                    closeFlag = 0;
                else
                    tolerance = obj.settings.fptolerance1;
                    closeFlag = 1;
                end
                if obj.config.guessdelta && ~isempty(obj.deltaJac)
                    newdelta = log(obj.edelta) + ...
                        obj.deltaJac*(rc_sigma - obj.oldsigma);
                    edel = obj.finddelta(rc_sigma, exp(newdelta), tolerance);
                else
                    edel = obj.finddelta(rc_sigma, obj.edelta, tolerance);
                end
                % Update oldsigma and edelta only in first stage and if
                % successful
                if closeFlag == 1 && max(isnan(edel)) < 1;
                   obj.oldsigma = rc_sigma;
                end   
                obj.edelta = edel;
                del = log(edel);
                
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
                            obj.deltaJac = obj.deltaJacobian(rc_sigma, obj.edelta);
                            g = 2* obj.deltaJac'* obj.ZWZ * xi;
                        end
                    else
                        xiX =  xi' * obj.X;
                        f = xiX * xiX';
                        if nargout > 1
                            obj.deltaJac = obj.deltaJacobian(rc_sigma, obj.edelta);
                            g = 2*obj.deltaJac'* obj.X * xiX';
                        end
                    end
                end
            end
            
            function stop = iterOutput(x, optimValues, state)
                stop = false;
                switch state
                    case 'iter'
                        fprintf('Iterations: %d / %d  \t Sigma: %2.4f\n', ...
                            obj.config.fp_invoc, obj.config.fp_iter, x(1));
                        obj.config.historyx = [obj.config.historyx ; x'];
                        obj.config.fp_iter = 0;
                        obj.config.fp_invoc = 0;
                end
            end
        end
        
        function R = estimate(obj, varargin)
            obj.edelta = exp(obj.ls);
            obj.initEstimation();
            obj.minimize([{'Display','iter-detailed'}, varargin{:}]);
            delta = log(obj.edelta);
            if isnan(delta)
                error('delta is NaN');
            end
            obj.rc_sigma = abs(obj.rc_sigma); % Works for symmetric or positive dist
            [obj.beta, obj.xi] = obj.lpart(delta);
            obj.alpha = obj.beta(strcmp(obj.getPriceName(), obj.vars));
            coef = [obj.beta; obj.rc_sigma];  % vector of all parameters
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
            obj.results.params.rc_sigma = obj.rc_sigma;
            R = obj.results.estimate;
        end
        
        function varcovar = computeVariance(obj)
            derdel = obj.deltaJacobian(obj.rc_sigma, obj.edelta);
            dgf = (size(obj.X, 1) - size(obj.X, 2));
            obj.results.params.dgf = dgf;
            if isempty(obj.Z)
                dertheta = [-obj.X derdel];
                varcovar = inv(dertheta' * dertheta);
                return
            end
            dertheta = [-obj.X derdel]'* obj.Z;
            if isfield(obj.config, 'test') && obj.config.test
                xiZ = bsxfun(@times, obj.xi, obj.Z );
                obj.W = inv(xiZ'*xiZ);
            end
            obj.results.other.cond = cond(dertheta * obj.W * dertheta');
            if obj.settings.optimalIV
                varcovar = inv(dertheta * obj.W * dertheta');
            else
                varcovar = (obj.xi' * obj.xi/dgf) .* inv(dertheta * obj.W * dertheta');
            end
        end
                      
%% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function init(obj, varargin)
            % init(data, selection) - both arguments optional
            if nargin > 1 && ~isempty(varargin{1})
                obj.T = varargin{1};
            end
            if nargin == 3 % To let initSimulation invoke init@DemandEstimate
                init@DemandEstimate(obj, varargin{2});
                return
            end

            init@DemandEstimate(obj);
            if nargin == 1
                obj.nonlinparams = [];
                obj.nonlintype = [];
                nonlin = {obj.var.nonlinear, obj.var.nonlinearlogs, ...
                    obj.var.nonlineartriangular};
                for i = 1:length(nonlin)
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
                if isempty(obj.rc_sigma)
                    if isempty(obj.config.randstream)
                        obj.rc_sigma = randn(length(obj.nonlinparams),1);
                    else
                        obj.rc_sigma = obj.config.randstream.randn(length(obj.nonlinparams),1);                    
                    end
                end
                if size(obj.rc_sigma, 2) > 1
                    obj.rc_sigma = obj.rc_sigma';
                end
                if length(obj.nonlinparams) ~= length(obj.rc_sigma)
                    error('rc_sigma and nonlinear have different dimensions');
                end
                obj.results.rc_sigma0 = obj.rc_sigma;
                obj.vars2 = ...
                    cellfun(@(x) {sprintf('rc_%s', x)}, obj.nonlinparams );
            end
            
            obj.x0 = obj.X(:,2:end); 
            obj.x2 = obj.T{:, obj.nonlinparams }; 
            nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
            if any(nonlinprice) % For CES, logprice rather than price.
                obj.x2(:, nonlinprice) = obj.getPrice();
            end
            obj.randdraws();
        end
        
        function randdraws(obj)
            K = size(obj.x2, 2); % Number of variables
            obj.v = [];
            obj.vx = [];
            if strcmpi(obj.settings.drawmethod, 'quadrature')
                obj.config.quadrature = true;
                [X, obj.quadw] = nwspgr('KPN', K, obj.settings.quaddraws);
                obj.settings.nind = length(obj.quadw);
                for k = 1:K
                    obj.v(:,:,k) = repmat(X(:,k)', size(obj.x2, 1), 1);
                    obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v(:,:,k));
                end
            else
                if isempty(obj.draws) % Can be set manually
                    if obj.settings.marketdraws
                        obj.draws = Sampling.draw(obj.settings.drawmethod, ...
                            K, obj.nmkt*obj.settings.nind, ...
                            obj.config.randstream);
                    else
                        obj.draws = Sampling.draw(obj.settings.drawmethod, ...
                            K, obj.settings.nind, obj.config.randstream);
                    end
                end
                for k = 1:K
                    % For each k, create draws for each market - wk (nind x nmkt)
                    if obj.settings.marketdraws
                        wk = reshape( obj.draws(:, k), obj.settings.nind, obj.nmkt )';
                    else
                        wk = repmat( obj.draws(:, k)', obj.nmkt, 1);
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
                    obj.vx(:,:,k) = bsxfun(@times, obj.x2(:,k), obj.v(:,:,k));
                end
            end
        end
        
        function initSimulation(obj, varargin)
            if nargin > 1
                selection = varargin{1};
                obj.init([], selection); 
                vnew = [];
                vxnew = [];
                for k = 1:size(obj.rc_sigma, 1)
                    vnew(:,:,k) = obj.v(selection,:,k);
                    vxnew(:,:,k) = obj.vx(selection,:,k);
                end
                obj.v = vnew;
                obj.vx = vxnew;
                obj.dummarket = dummyvar(obj.marketid);
            end
            obj.x0 = obj.X(:,2:end); % Remove prices from X
            obj.x2 = obj.T{:, obj.nonlinparams };
            nonlinprice = strcmp(obj.var.price, obj.nonlinparams);
            if any(nonlinprice) % For CES, logprice rather than price.
                obj.x2(:, nonlinprice) = obj.getPrice();
            end
            obj.nlpart(obj.rc_sigma);
            if ~isempty(obj.beta) && ~isempty(obj.rc_sigma)
                % Create starting values for finddelta
                if isempty(obj.s)
                    obj.s = sum(bsxfun(@times, obj.dummarket, ...
                        1 ./ sum(obj.dummarket) ), 2) .* 0.5;
                    ed = obj.s .* 2;
                else
                    ed = obj.s ./ obj.s0;
                end
                obj.edelta = obj.finddelta(obj.rc_sigma, ed, obj.settings.fptolerance1);
                obj.alpha = obj.beta(1);
                obj.nlpart(obj.rc_sigma)
                obj.d = log(obj.edelta) - obj.alpha*obj.getPrice(); 
            else
                error('MixedLogitDemand.initSimulation needs beta and rc_sigma to be set');
            end
        end
            
        function createPeriods(obj)
            obj.period = cell(obj.nmkt,1);
            periodCopy = copy(obj); 
            periodCopy.T = [];
            for t = 1:max(obj.marketid)
                selection = obj.marketid == t;
                newperiod = copy(periodCopy); 
                obj.periodSelect(newperiod, selection);
                obj.period{t} = newperiod;
            end            
        end
                
        function periodSelect(obj, newperiod, selection)
            newperiod.X = obj.X(selection, :);
            newperiod.p = obj.p(selection, :);
            newperiod.ms = obj.ms(selection, :);
            newperiod.s = obj.s(selection, :);
            newperiod.marketid = obj.marketid(selection, :);
            newperiod.T = obj.T(selection, :);
            
            newperiod.v = obj.v(selection, :, :);
            newperiod.vx = obj.vx(selection, :, :);
            newperiod.x0 = obj.x0(selection, :);
            newperiod.x2 = obj.x2(selection, :);
            newperiod.edelta = obj.edelta(selection, :);
            newperiod.d = obj.d(selection, :);
            
            if ~isempty(obj.dummyvars)
                newperiod.dummyvars = obj.dummyvars(selection, :);
            end
            newperiod.dummarket = dummyvar(obj.marketid(selection, :));
        end
                
        function demand = pack(obj, filename)
            demand = MixedLogitDemand;
            demand.settings = obj.settings;
            demand.var = obj.var;
            demand.config = obj.config;
            demand.results = obj.results;
            demand.beta = obj.beta;
            demand.rc_sigma = obj.rc_sigma;
        end
        
        function obj = MixedLogitDemand(varargin)
            if nargin > 0
                obj.T = varargin{1};
            end
            obj.var.nonlinear = ''; % non-linear parameters other than price
            obj.var.nonlinearlogs = ''; % log-normal parameters
            obj.var.nonlineartriangular = ''; % triangular dist parameters

            obj.settings.optimalIV = false;
            obj.settings.drawmethod = 'hypercube';
            obj.settings.marketdraws = false; % Different draws for each market
            obj.settings.nind = 100; % number of simulated "indviduals" per market 
            obj.settings.quaddraws = 10;
            obj.settings.maxiter = 100;
            obj.settings.fptolerance1 = 1e-14; % lower tol for first FP its
            obj.settings.fptolerance2 = 1e-14; % high tol for last iterations
            obj.settings.ces = false;
            obj.config.tolerance = 1e-6;
            obj.config.fpmaxit = 1000; % maximum iterations in contraction mapping
            obj.config.restartMaxIterations = 5; % Max # of restarts if fval is high
            obj.config.restartFval =  10^3; % Restart optimal IV estimation if fval >.
            obj.config.quadrature = false; % Internal setting
            obj.config.analyticJacobian = true;
            obj.config.randstream = []; % random stream for estimation in parallel 
            obj.config.estimateDescription = 'Random Coefficient Logit Demand'; 
            obj.config.guessdelta = true;
            obj.config.fp_iter = 0;
            obj.config.fp_invoc = 0;
        end      
    end

%% Basic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    methods % (Access = private)
        function newedelta = finddelta(obj, rc_sigma, edelta, tolerance)
            selection = logical(obj.dummarket);
            obj.nlpart(rc_sigma); % Creates obj.mu
            newedelta = edelta;
            for t = 1:max(obj.marketid)
                i = 0;
                maxdev = 100;
                edeltasel = edelta(selection(:,t),1);
                expmusel = obj.expmu(selection(:,t), :);
                sorig = obj.s(selection(:,t),1);
                while maxdev > tolerance && i < obj.config.fpmaxit
                    eg = bsxfun(@times, edeltasel, expmusel );
                    denom = 1 ./ ( 1 + sum(eg));
                    if obj.config.quadrature
                        share = bsxfun(@times, eg, denom) * obj.quadw; % Row means
                    else
                        share = mean(bsxfun(@times, eg, denom), 2); % Row means
                    end
                    newedel = edeltasel .* sorig ./ share;
                    maxdev = max(abs(newedel - edeltasel));
                    edeltasel = newedel;
                    i = i + 1;
                end
                newedelta(selection(:,t)) = edeltasel;
            end
        end
        
        function [ bet, xi ] = lpart(obj, del )
            if ~isempty(obj.Z)
                bet = obj.inv_x1ZWZx1 * (obj.X' * obj.ZWZ * del);
            else
                bet = obj.inv_x1ZWZx1 * (obj.X' * del);
            end
            xi = del - obj.X * bet;
        end
        
        function nlpart(obj, rc_sigma)
        % Computes the mu or nonlinear part of utility
            [n,m,K] = size(obj.vx);
            mu = zeros(n, m);
            for k = 1:K
                mu = mu + obj.vx(:,:,k) * rc_sigma(k);
            end
            obj.expmu = exp(mu);
        end
        
        function [f,indshare] = share(obj, edelta)
            eg = bsxfun(@times, edelta, obj.expmu );
            denom1 = 1 ./ ( 1 + obj.dummarket' * eg);
            denom = denom1(obj.marketid,:);
            indshare = eg .* denom;
            if obj.config.quadrature
                f = indshare * obj.quadw; 
            else
                f = mean(indshare, 2); 
            end
        end

        function initEstimation(obj)
            if ~obj.settings.optimalIV
                if isempty(obj.W) && ~isempty(obj.Z) % W can be specified manually
                    obj.W = inv(obj.Z' * obj.Z);
               end
            else
                if ~isempty(obj.optimalInstrument )
                    return % Already calculated optimal IV
                end
                pHat = obj.Z*((obj.Z'*obj.Z)\obj.Z'*obj.getPrice()); 
                deltaHat = [pHat, obj.x0] * obj.beta;
                obj.optimalInstrument = obj.deltaJacobian(obj.rc_sigma, ...
                    exp(deltaHat));
                obj.Z  = [obj.x0 pHat obj.optimalInstrument];
                xiZ = bsxfun(@times, obj.xi, obj.Z );
                obj.W = inv(xiZ'*xiZ);
            end
            
            if ~isempty(obj.Z)
                obj.ZWZ = obj.Z * obj.W * obj.Z';
                obj.inv_x1ZWZx1 = inv(obj.X'*obj.ZWZ*obj.X);
            else
                obj.inv_x1ZWZx1 = inv(obj.X'*obj.X);
            end
        end
    end
end

