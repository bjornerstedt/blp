classdef RCDemand < NLDemand
    % Random Coefficient demand class
    %   Simulation of shares and estimation of parameters beta and sigma.
   
     properties (SetAccess = protected, Hidden = true )
        draws % Random draws (nind, nmkt, k) matrix k = # random params. 
        W
        x2 % nonlinear parameter vector
        xi % Unobservable utility, residual
        vars2 % Nonlinear variable names in output (with rc_ prefix)
        
        edelta % Saved between invocations of objective
        deltaJac % Used to guess new edelta
        int % Internal variables
    end
    
    methods
%% GENERAL DEMAND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function [s, si] = shares(obj, p, market)
            [s, si] = obj.period{market}.shares( p );
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
       
        function sh = shareJacobian(obj, P, market)
            sh = obj.period{market}.shareJacobian(P);          
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
            obj.int.oldsigma = zeros(size(obj.sigma));
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
            if fval == 1e10
                error('RCDemand.estimate could not find an estimate')
            end
            obj.results.other.fval = fval;
            obj.results.other.restarts = i - 1;
            obj.results.other.minimum = finished;
        end
        
        function [f, g] = objective(obj, sigma)
            % This function defines the objective over which to minimize
            obj.edelta =  obj.findDelta(sigma);
            del = log(obj.edelta);
            
            if ~isreal(del) || max(isnan(del)) == 1
                f = 1e+10;
                if nargout > 1
                    g = 1e+10 * ones(size(obj.sigma));
                end
                disp('Negative deltas found at observations')
                err = find(isnan(del));
                if length(err) < 11
                    disp(err)
                end
            else
%             if strcmpi(obj.settings.paneltype, 'fe')
%                 davt = accumarray(obj.panelid, del, [], @mean);
%                 davt = davt(obj.panelid, :);
%                 del = del - davt;
%             end
                xi = obj.int.annihalator * del;
                if ~isempty(obj.Z)
                    f = xi' * obj.int.ZWZ * xi;
                    if nargout > 1
                        obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
                        g = 2 * obj.deltaJac'* obj.int.ZWZ * xi;
                    end
                else
                    f = xi' * xi;
                    if nargout > 1
                        obj.deltaJac = obj.deltaJacobian(sigma, obj.edelta);
                        g = 2 * obj.deltaJac' * xi;
                    end
                end
            end
        end
        
        function truevals = setTrueResults(obj, beta)
            beta = [-obj.alpha; beta'];
            if strcmpi(obj.settings.paneltype, 'fe')
                beta = beta(1:end-1);
            end
            truevals = table([beta; obj.sigma]);
            truevals.Properties.VariableNames = {'True_val'};
            obj.results.trueValues = truevals;
        end
        
        function R = estimate(obj, varargin)
            obj.init(varargin{:});
            obj.edelta = exp(obj.share.ls);
            if obj.settings.optimalIV
                if isempty(obj.Z)
                    error('RCDemand.settings.ptimalIV only works if instruments are set.')
                end
                est = obj.estimationStep(false, varargin{:});
                R = obj.estimationStep(true, varargin{:});
                obj.results.estimate1 = est;
            else
                R = obj.estimationStep(false, varargin{:});
            end
            % Create starting values for findDelta
            obj.d = log(obj.edelta) + obj.alpha * obj.Xorig(:, 1);
            for t = 1:max(obj.marketid)
                obj.period{t}.d = obj.d(obj.dummarket(:, t));
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
            obj.beta = obj.int.estimationMatrix * delta;
            xi = delta - obj.X * obj.beta;
            
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
            obj.results.alpha = obj.alpha;
            obj.results.sigma = obj.sigma';
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
            if isfield(obj.results,'trueValues') 
                obj.results.estimate = [obj.results.trueValues, obj.results.estimate];
            end        
            R = obj.results.estimate;
        end
        
        function varcovar = computeVariance(obj)
            derdel = obj.deltaJacobian(obj.sigma, obj.edelta);
            if isempty(obj.Z)
                dertheta = [-obj.X derdel];
                invXX = inv(dertheta' * dertheta);
                if obj.settings.robust
                    xiX =bsxfun(@times, obj.xi, dertheta);
                    varcovar = invXX * (xiX' * xiX) * invXX;
                    % TODO: Add dgf adjustment
                else
                    varcovar = ((obj.xi' * obj.xi) ./ obj.results.dgf) * invXX;
                end
            else
%                 if obj.config.test
%                     pHat = obj.Z*((obj.Z'*obj.Z)\obj.Z'*obj.X(:, 1));
%                     deltaHat = [pHat, obj.X(:,2:end)] * obj.beta;
%                 end

                dertheta = [obj.X, derdel]'* obj.Z;
% Robust estimation if test is set:
%             if isfield(obj.config, 'test') && obj.config.test
%                 xiZ = bsxfun(@times, obj.xi, obj.Z );
%                 obj.W = inv(xiZ'*xiZ);
%             end
                obj.results.other.cond = cond(dertheta * obj.W * dertheta');
                if obj.settings.robust
                    varcovar = inv(dertheta * obj.W * dertheta');
                else
                    varcovar = (obj.xi' * obj.xi / obj.results.dgf) .* ...
                        inv(dertheta * obj.W * dertheta');
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
                nonlinparams = obj.randdraws();
            else
                nonlinparams = obj.draws.nonlinparams;
            end
            obj.x2 = obj.data{:, nonlinparams };
            % log price in CES in x2:
            nonlinprice = strcmp(obj.var.price, nonlinparams);
            if any(nonlinprice) && obj.settings.ces
                obj.x2(:, nonlinprice) = log(obj.data{:, obj.var.price});
                nonlinparamsCES = nonlinparams;
                nonlinparamsCES{nonlinprice} = 'lP';
                obj.vars2 = ...
                    cellfun(@(x) {sprintf('rc_%s', x)}, nonlinparamsCES );
            else
                obj.vars2 = ...
                    cellfun(@(x) {sprintf('rc_%s', x)}, nonlinparams );
            end
            if ~isempty(selection)
                obj.x2 = obj.x2(selection, :); 
            end
            obj.getSigma();
            obj.initPeriods();
        end
        
        function nonlinparams = randdraws(obj)
            if obj.settings.marketDraws
                markets = max(obj.marketid);
            else
                markets = 1;
            end
            obj.draws = Draws('DrawMethod', obj.settings.drawMethod,...
                'Markets', markets, 'Individuals', obj.settings.nind,...
                'Accuracy', obj.settings.accuracy, ...
                'RandStream', obj.config.randstream);
            if isempty(obj.var.nonlinear)
                error('Some variable has to be specified as nonlinear');
            end
            nonlinparams = obj.draws.create( obj.var.nonlinear);
        end
        
        function getSigma(obj)
        % obj.sigma is set in simulation, in estimation it is created
        % from obj.settings.sigma0 or from random draw. 
            nonlinparams = obj.draws.nonlinparams;
            if ~isempty(obj.settings.sigma0)
                obj.sigma = obj.settings.sigma0;
                obj.results.sigma0 = obj.sigma;
                if size(obj.sigma, 2) > 1
                    obj.sigma = obj.sigma';
                end
                if length(nonlinparams) ~= length(obj.sigma)
                    error('sigma and nonlinear have different dimensions');
                end
            end
            if isempty(obj.sigma)
                if isempty(obj.config.randstream)
                    obj.sigma = randn(length(nonlinparams), 1);
                else
                    obj.sigma = obj.config.randstream.randn(length(nonlinparams),1);
                end
                obj.results.sigma0 = obj.sigma;
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
                if obj.settings.marketDraws
                    newmarket.v = obj.draws.draws(:,:,t)'; 
                else
                    newmarket.v = obj.draws.draws'; 
                end
                newmarket.init();
                obj.period{t} = newmarket;
            end
        end
              
        function obj = RCDemand(varargin)
            obj = obj@NLDemand(varargin{:});
            obj.var.setParameters({'nonlinear'});
            obj.settings.setParameters({'optimalIV','drawMethod',... 
                'marketDraws','nind','accuracy','maxiter', 'sigma0'});
            obj.config = SettingsClass({'tolerance','fptolerance1','fptolerance2', ...
                'restartMaxIterations','restartFval', 'test', 'fpmaxit',...
                'randstream','hessian','guessDelta','compiled','quietly'});
            
            obj.settings.paneltype = 'lsdv';

            obj.settings.optimalIV = false;
            obj.settings.drawMethod = 'hypercube';
            obj.settings.marketDraws = false; % Different draws for each market
            obj.settings.nind = 100; % number of simulated "indviduals" per market 
            obj.settings.accuracy = 10;
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
            obj.config.guessDelta = true;
            obj.config.hessian = false;
            obj.config.compiled = false; % Use c++ code
            obj.config.quietly = true;
            obj.results.estimateDescription = 'Random Coefficient Logit Demand'; 
            obj.results.sigma0 = []; 
        end  
        
    end

%% Basic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    methods % (Access = private)
        function newedelta = findDelta(obj, sigma)
            if max(abs(sigma - obj.int.oldsigma)) < 0.01
                tolerance = obj.config.fptolerance2;
                closeFlag = 0;
            else
                tolerance = obj.config.fptolerance1;
                closeFlag = 1;
            end
            if obj.config.guessDelta && ~isempty(obj.deltaJac)
                newdelta = log(obj.edelta) + ...
                    obj.deltaJac * (sigma - obj.int.oldsigma); 
                edelta = exp(newdelta);
            else
                edelta =  obj.edelta;
            end
            sel = logical(obj.dummarket);
            newedelta = zeros(size(edelta));
            nccount = 0;
            for t = 1:max(obj.marketid)
                [newedelta(sel(:,t)), nonconv] = obj.period{t}.findDelta(sigma, ...
                    edelta(sel(:,t)), tolerance);
                nccount = nccount + nonconv;
            end
            if nccount > 0
                warning('findDelta did not converge in %d markets', nccount)
            end
           % Update oldsigma and edelta only in first stage and if
            % successful
            if closeFlag == 1 && max(isnan(newedelta)) < 1;
                obj.int.oldsigma = sigma;
            end
        end       
            
        function initEstimation(obj, optIV)
        % initEstimation is invoked by estimateStep
            if ~optIV
                if isempty(obj.W) && ~isempty(obj.Z) % W can be specified manually
                    obj.W = inv(obj.Z' * obj.Z);
                end
                obj.int.isOptimalIV = false; 
            else
                if obj.int.isOptimalIV
                    return % Already calculated optimal IV
                end
                obj.int.isOptimalIV = true;
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
                obj.W = inv(xiZ' * xiZ);
            end
            
            if ~isempty(obj.Z)
                obj.int.ZWZ = obj.Z * obj.W * obj.Z';
                obj.int.estimationMatrix = inv(obj.X'*obj.int.ZWZ*obj.X) * obj.X' * obj.int.ZWZ;
            else
                obj.int.estimationMatrix = inv(obj.X'*obj.X) * obj.X';
            end
            obj.int.annihalator = eye(size(obj.X, 1)) - obj.X * obj.int.estimationMatrix;
        end
        
    end
end

