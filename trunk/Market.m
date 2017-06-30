classdef Market < Estimate 
    % MARKET Calculates market equilibrium
    %   Based on ownership structure and conduct.
    
    properties
        firm
        p %prices
        c %Costs
        gamma = 0 % Scale effects
        demand % Demand class        
        marketid % Protected?
    end
    properties (SetAccess = protected, Hidden = true )
        RR %ownership tranformation with conduct
    end
    
    methods
       function init(obj)
            if isempty(obj.demand)
                error('Demand object must be specified')
            end
            % Removes selection from demand:
            obj.demand.init();
            if isempty(obj.var.firm)
                if isempty(obj.firm)
                    error('var.firm has to be specified');
                end
                % obj.firm can be set directly if obj.var.firm is empty.
            else
                obj.firm = obj.demand.data{:, obj.var.firm};
                if ~iscategorical(obj.firm)
                    [~, ~, obj.firm] = unique(obj.firm);
                end
            end
            if ~isempty(obj.data)
                init@Estimate(obj);
                obj.panelid = obj.demand.panelid;
            end
            % Use weighting of averages as in demand unless set by user
            if isempty(obj.settings.valueShares)
                obj.settings.valueShares = obj.demand.useValueShares();
            end
            % Copied, can be null
            obj.p = obj.demand.data{:, obj.demand.var.price};
            obj.marketid = obj.demand.marketid;
        end
        
        function initSimulation(obj, market_number)
            R = dummyvar(obj.firm(obj.marketid == market_number) )';
            obj.RR  = R' * R;
            if obj.settings.conduct > 0
                obj.RR(obj.RR == 0) = obj.settings.conduct;
            end
        end
        
        function exitflag = equilibrium(obj, varargin)
        % equilibrium([selection]) calculates the equilibrium price for selection
        % The new price is put in the Market.p property. 
            obj.init();
            if nargin > 1
                selection = varargin{1} & ~isnan(obj.c);
            else
                selection = ~isnan(obj.c);
            end           
            if ~isempty(obj.demand.share)
                p0 = obj.p;
            else
                p0 = obj.c;                
            end
            marketid_list = unique(obj.marketid(selection));
            options = optimoptions(@fsolve,... 
                'MaxFunEvals',obj.settings.maxit,'Display', 'off');
            obj.p = nan(size(obj.marketid));
            convCount = 0;
            for i = 1:length(marketid_list)
                t = marketid_list(i);
                selection = obj.marketid == t;
                obj.initSimulation(t);
                [pt, ~, exitflag, output] = fsolve( ...
                    @(x)obj.foc(x, obj.c(selection), t), ...
                    p0(selection), options);
                obj.p(selection) = pt;
                obj.results.iterations = output.iterations;
                if exitflag == 1
                    convCount = convCount+1;
                end
            end
            if ~obj.config.quietly
                display(sprintf('Simulation converged for %d of %d markets', ...
                    convCount, length(marketid_list)));
            end
            obj.results.convCount = convCount;
        end
  
        function findCosts(obj, varargin)
        % findCosts([selection]) calculates costs for selection
            obj.init();
            if nargin > 1
                selection = varargin{1};
                marketid_list = unique(obj.marketid(selection));
            else
                selection = ones(size(obj.marketid));
                marketid_list = unique(obj.marketid);
            end
            quantity = obj.demand.actualDemand();
            obj.c = nan(size(obj.marketid));
            obj.results.findCosts.cond = [];
            for i = 1:length(marketid_list)
                t = marketid_list(i);
                msel = obj.marketid == t;
                obj.initSimulation(t);
                [sj, cnd] = linsolve( obj.RR .* obj.demand.shareJacobian(obj.p(msel), t), ...
                    -quantity(msel) );
                obj.c(msel) = obj.p(msel) - sj;
                obj.results.findCosts.cond = min([obj.results.findCosts.cond, cnd]);
            end
            obj.results.findCosts.negCostCount = sum(obj.c(selection) < 0);
            if obj.results.findCosts.negCostCount < 0
                warning('Calculated negative costs')
            end
            if min(obj.p(selection) - obj.c(selection)) < 0
                warning('Calculated costs exceeding prices')
            end
        end

		function f = foc(obj,  P, ct, t)
			S = obj.demand.shares( P, t );
            ct = ct + obj.gamma * S; % Scale effects
			f = ( obj.RR .* obj.demand.shareJacobian( P, t )) * (P - ct) + S;
		end 

        function theta = estimate(obj, varargin)
            obj.data = obj.demand.data;
            obj.init(varargin{:});
            theta = estimate@Estimate(obj, varargin{:});
 		end 

        function R = estimateGMM(obj, theta)
            % Simultaneous 2SLS estimate of demand and costs over alpha
            % obj.demand.residuals() calculates residuals after OLS,
            % taking alpha and sigmas as gvien in theta, along the lines of
            % BLP estimation.
            obj.data = obj.demand.data;
            obj.init();
            obj.demand.estimate();
            obj.findCosts();
            obj.y = obj.c;
            obj.estimate();
            options = optimoptions(@fminunc, 'Display', 'iter' ,'Algorithm', 'quasi-newton', 'MaxIter',50);
            
            Z = [bsxfun(@times, obj.demand.Z, obj.demand.residuals(theta)), ...
                bsxfun(@times, obj.X, obj.results.xi)];
            W = inv(Z' * Z);
% Test with quadratic optimization
%             C = chol(W);
            [theta] = fminunc(@(x)objective(x, obj), theta, options);
%             [theta] = lsqnonlin(@(x)objective(x, obj), theta);
            
            Z = [bsxfun(@times, obj.demand.Z, obj.demand.results.xi), ...
                bsxfun(@times, obj.X, obj.residuals())];
            W = inv(Z' * Z);
            [theta] = fminunc(@(x)objective(x, obj), theta, options);
%             [theta] = lsqnonlin(@(x)objective(x, obj), theta);

            [~, beta ] = obj.demand.residuals(theta);
            [~, lambda ] = obj.residuals();
            theta = [beta; lambda];
            Sxz = [obj.demand.X, zeros(size(obj.X)); ...
                zeros(size(obj.demand.X)), obj.X ]' * ...
                [obj.demand.Z, zeros(size(obj.X)); ...
                zeros(size(obj.demand.Z)), obj.X ];
            varcovar = inv(Sxz * W * Sxz');
            se = sqrt(  diag(varcovar) );  
            tvalue = theta ./ se;
            R = table(theta,se,tvalue);
            function val = objective(theta, obj)
                % Residuals as function of demand params:
                xi = obj.demand.residuals(theta);
                obj.findCosts();
                eta = obj.residuals();
                
                xiZ = xi' * obj.demand.Z;
                etaZ = eta' * obj.X;
                val = [xiZ, etaZ] * W * [xiZ, etaZ]';
%                 val = C * [xiZ, etaZ]';
            end
        end
        
        function [xi, gamma] = residuals(obj)
            gamma = (obj.X' * obj.X) \ ( obj.X' * obj.c);
            xi = obj.c - obj.X * gamma;
        end
        
		function R = estimateCosts(obj)
            obj.init();
            obj.findCosts();
            if min(obj.c) > 0
                lc = log( obj.c );
            else
                warning('Calculated negative costs');
                lc = log( max(obj.c, .0001) );
            end
            obj.y = lc;
            R = obj.estimate();
        end
        
        function res = getMarketShares(obj)
            % getMarketShares([weights]) Market shares by market.
            % Used as weights in summary and compare
            % When should actual and simulated prices be used? 
            q = obj.demand.getDemand(obj.p);
            if obj.settings.valueShares % Note that valueShares should be true for ces.
                qu = q * obj.p;
            else
                qu = q;
            end
            qbar = accumarray(obj.marketid, qu, [], @sum);
            qbar = qbar(obj.marketid, :);
            res = qu ./ qbar;
        end
        
        function markettab = summary(obj, varargin)
            % SUMMARY(options) creates summary by varargin variable
            % Weighted or unweighted means over markets of: p, c and Lerner index.
            % Within group, means can be unweighted or weighted by quantity
            % or by value. Without group_vars, reporting is at the product level.
            % Output is as a table.
            args = inputParser;
%             args.addOptional('obj2', []);
            args.addOptional('obj2', [], @(x)isa(x, 'Market'));
            args.addParameter('selection',[], @islogical);
            args.addParameter('GroupingVariables', 'Firm', @ischar);
            args.addParameter('weights', [], @ischar);
            args.addParameter('WeightedAverages', ... 
                obj.settings.weightedAverages, @islogical);
            args.parse(varargin{:});
            if ~strcmpi(args.Results.GroupingVariables, 'Firm')
                indexvars = {args.Results.GroupingVariables, obj.var.market};
                var = obj.demand.data(:, indexvars);
            else
                indexvars = {'Firm', obj.demand.var.market};
                var = table(obj.firm, obj.marketid);
                var.Properties.VariableNames = indexvars;
            end
            
            if isempty(args.Results.obj2)
                %%%%%% Single market summary
                if isempty(args.Results.selection)
                    selection = ~isnan(obj.c);
                else
                    selection = args.Results.selection & ~isnan(obj.c);
                end
                res = [var, array2table([obj.p, obj.c]) ];
                res.Lerner = (obj.p - obj.c) ./ obj.p ;
                tableCols = {'Price', 'Costs', 'Lerner'};
                res.Properties.VariableNames = [indexvars, tableCols];
                shares = true;
            else
                %%%%%% Two market comparison
                obj2 = args.Results.obj2;
                if isempty(args.Results.selection)
                    selection = ~isnan(obj2.p);
                else
                    selection = args.Results.selection & ~isnan(obj2.p);
                end
                tableCols = {'Costs', 'Price1', 'Price2', 'PriceCh'};
                priceChange = (obj2.p - obj.p ) ./ obj.p;
                res = [var, array2table([ obj.c, obj.p,  obj2.p, priceChange],...
                    'VariableNames', tableCols)];
                shares = false;
            end
            
            % weights can be set manually
            if ~args.Results.WeightedAverages
                weights = [];
            elseif ~isempty(args.Results.weights)
                weights = args.Results.weights;
            else
                weights = obj.getMarketShares();
            end
            if ~all(selection)
                res = res(selection, :);
                if ~isempty(weights)
                    weights = weights(selection, :);
                end
            end
            x = Estimate.means(res, tableCols, indexvars, weights);
            if shares
                sh = obj.getMarketShares();
                [~, ~, rowIdx] = unique(res(:, indexvars), 'rows');
                x.MarketSh = accumarray(rowIdx, sh(selection));
                tableCols = [tableCols, {'MarketSh'}];
                markettab = Estimate.means(x, tableCols, args.Results.GroupingVariables);
            else
                markettab = Estimate.means(x, tableCols, args.Results.GroupingVariables);
            end
        end
                       
        function obj = Market(varargin)
            obj = obj@Estimate();
            if nargin > 0
                obj.demand = copy(varargin{1});
            end
            varsEstimate = {'market','panel','depvar','exog', ...
                'endog','instruments'};
            obj.var = SettingsClass([varsEstimate, {'firm'}]);
            obj.settings.setParameters({'conduct', 'dampen', 'maxit', ...
                'weightedAverages', 'valueShares'});
            obj.settings.conduct = 0;  
            obj.settings.maxit = 3000;
            obj.settings.dampen = 1;
            obj.settings.weightedAverages = true;
            obj.settings.valueShares = false;
        end      
        
		function f = margins(obj,  P)
			S = obj.demand.shares( P );
			f = linsolve( (-obj.RR) .* (obj.demand.shareJacobian( S , P)) , S );
        end
        
        function f = fixedPoint(obj, maxit)
			convergence = 1 ;
            if ~isempty(obj.demand.share)
                P = obj.p;
            else
                P = obj.c;                
            end
            sensitivity = 10^-6;
			dist = sensitivity + 1;
			i = 0;
			diff = 0;
			for i = 1:maxit 
                if dist < sensitivity
                    break;
                end
				%for linear logit set Q=S, for CES logit set Q=S./P xxx
				Pn = obj.c + obj.margins( P);
				if isnan(max(Pn))
					i = maxit;
					convergence = 0;
					warning( 'Could not invert the share Jacobian.' );
					break;
				end 
				diff = Pn - P;
				dist = dot(diff, diff);
		% 	diff = foc(obj.D, P);
		% 	dist = dot(diff, diff);
				if min(P) < 0 
					convergence = 0;
					warning( 'Negative prices in price vector' );
					break;
				end 			
				P = obj.settings.dampen * Pn + (1 - obj.settings.dampen) * P;
			end 
			if i >= maxit 
				warning( 'Max number of iterations exceeded.' );
				convergence = 0;
			end 
			S = obj.demand.shares( P ); %This function should be called demand
			diff = dot(S, S);
			if diff < sensitivity 
				warning('Converging to small shares.' );
				convergence = 0;
			end 
			if min(P) < 0 || min(S) < 0  
				warning( 'Negative prices or shares.' );
				convergence = 0;
			end 
			%st_numscalar('r(maxpricediff)', max(abs(diff)) );
			diff = obj.foc( P);
			%st_numscalar('r(fixedpointdiff)', cross(diff, diff) )	;
			obj.p = P;
			f = [i,convergence];
        end
        
        function newmarket = clone(obj)
            newmarket = copy(obj);
        end      
    end
end

