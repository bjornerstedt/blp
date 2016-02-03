classdef Market < Estimate 
    % MARKET Calculates market equilibrium
    %   Based on ownership structure and conduct.
    
    properties
        firm
        q %quantities
        p %prices
        p0 %Initial prices
        c %Costs
        gamma = 0 % Scale effects
        marketid % Protected?
        demand % Demand class        
    end
    properties (SetAccess = protected, Hidden = true )
        s %shares
        RR %ownership tranformation with conduct
        c1, p1, q1 %Views for bootstrapping
        sel % market class with selection
        sim % per market data, set in initSimulation
    end
    
    methods
       function init(obj)
            if isempty(obj.demand)
                error('Demand object must be specified')
            end
            obj.demand.init();
            if isempty(obj.var.firm) 
                error('var.firm has to be specified');
            end
            % Set unless new firm has been set for example as post merger
            % ownership.
            if isempty(obj.firm)
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
            obj.p = obj.demand.p;
            obj.q = obj.demand.q; 
            obj.marketid = obj.demand.marketid;
            % With simulated demand, obj.share has not been set
            if ~isempty(obj.demand.share)
                obj.p0 = obj.p;
                obj.s = obj.demand.actualDemand();
            else
                obj.p0 = obj.c;                
            end
        end
        
        function initSimulation(obj, market_number)
            obj.demand.initSimulation(market_number);
            R = dummyvar(obj.firm(obj.marketid == market_number) )';
            obj.RR  = R' * R;
            if obj.settings.conduct > 0
                obj.RR(obj.RR == 0) = obj.settings.conduct;
            end
        end
        
        function exitflag = equilibrium(obj, varargin)
            obj.init();
            if nargin > 1
                selection = varargin{1} & ~isnan(obj.c);
            else
                selection = ~isnan(obj.c);
            end           
            marketid_list = unique(obj.marketid(selection));
            options = optimoptions(@fsolve,... 
                'MaxFunEvals',obj.settings.maxit,'Display', 'off');
            obj.p = nan(size(obj.marketid));
            obj.s = nan(size(obj.marketid));
            convCount = 0;
            for i = 1:length(marketid_list)
                t = marketid_list(i);
                selection = obj.marketid == t;
                obj.initSimulation(t);
                [pt, ~, exitflag, output] = fsolve( ...
                    @(x)obj.foc(x, obj.c(selection)), ...
                    obj.p0(selection), options);
                obj.p(selection) = pt;
                obj.s(selection) = obj.demand.shares(pt);
                obj.results.iterations = output.iterations;
                if exitflag == 1
                    convCount = convCount+1;
                end
            end
            display(sprintf('Simulation converged for %d of %d markets', ...
                convCount, length(marketid_list)));
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
            obj.c = nan(size(obj.marketid));
            obj.results.findCosts.cond = [];
            for i = 1:length(marketid_list)
                t = marketid_list(i);
                msel = obj.marketid == t;
                obj.initSimulation(t);
                [sj, cnd] = linsolve( obj.RR .* obj.demand.shareJacobian([]), ...
                    -obj.s(msel) );
                obj.c(msel) = obj.p(msel) - sj;
                obj.results.findCosts.cond = min([obj.results.findCosts.cond, cnd]);
            end
            obj.results.findCosts.minCosts = sum(obj.c(selection) < 0);
            if obj.results.findCosts.minCosts < 0
                warning('Calculated negative costs')
            end
            if min(obj.p(selection) - obj.c(selection)) < 0
                warning('Calculated costs exceeding prices')
            end
        end

		function f = foc(obj,  P, ct)
			S = obj.demand.shares( P );
            ct = ct + obj.gamma * S; % Scale effects
			f = ( obj.RR .* obj.demand.shareJacobian( P )) * (P - ct) + S;
		end 

        function theta = estimate(obj, varargin)
            obj.data = obj.demand.data;
            obj.init(varargin{:});
            theta = estimate@Estimate(obj, varargin{:});
 		end 

        function theta = estimateGMM(obj, theta)
            % Simultaneous 2SLS estimate of demand and costs over alpha
            obj.data = obj.demand.data;
            obj.init();
            WC = inv(obj.X' * obj.X);
            W = inv(obj.demand.Z' * obj.demand.Z);
            options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'MaxIter',50);
            
            [theta] = fminunc(@(x)objective(x, obj), theta, options);
            [~, beta ] = obj.demand.residuals(theta);
            [~, lambda ] = obj.residuals();
            theta = [beta; lambda];
            
            function val = objective(theta, obj)
                % Residuals as function of demand params:
                xi = obj.demand.residuals(theta);
                obj.findCosts();
                eta = obj.residuals();
                
                xiZ = xi' * obj.demand.Z;
                etaZ = eta' * obj.X;
                val = xiZ * W * xiZ' + etaZ * WC * etaZ';
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
        
        function res = getMarketShares(obj, varargin)
            if nargin > 1 && ~isempty(varargin{1})
                qu = varargin{1};
            elseif obj.settings.valueShares % Note that valueShares should be true for ces.
                qu = obj.q * obj.p;
            else
                qu = obj.q;
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
            args.addParameter('selection',[], @islogical);
            args.addParameter('GroupingVariables', 'Firm', @ischar);
            args.addParameter('weights', [], @ischar);
            args.parse(varargin{:});
            if ~strcmpi(args.Results.GroupingVariables, 'Firm')
                varnames = {args.Results.GroupingVariables, obj.var.market};
                var = obj.demand.data(:, varnames);
            else
                varnames = {'Firm', obj.demand.var.market};
                var = table(obj.firm, obj.marketid);
                var.Properties.VariableNames = varnames;
            end
            
            %%%%%% Specific code
            res = [var, array2table([obj.p, obj.c]) ];
            res.Lerner = (obj.p - obj.c) ./ obj.p ;
            tableCols = {'Price', 'Costs', 'Lerner'};
            
            if isempty(args.Results.selection)
                selection = ~isnan(obj.c);
            else
                selection = args.Results.selection & ~isnan(obj.c);
            end
            res.Properties.VariableNames = [varnames, tableCols];
            %%%%%% End specific code
            
            % weights can be set manually
            if ~obj.settings.weightedAverages
                weights = [];
            else
                weights = obj.getMarketShares(args.Results.weights);
            end
            if ~all(selection)
                res = res(selection, :);
                if ~isempty(weights)
                    weights = weights(selection, :);
                end
            end
            x = Estimate.means(res, tableCols, varnames, weights);
            markettab = Estimate.means(x, tableCols, args.Results.GroupingVariables);
        end
                     
        function [pricetab] = compare(obj, obj2, varargin)
            args = inputParser;
  %          args.addRequired('obj2', @(x)isa(x, 'NestedLogitDemand'));
            args.addParameter('selection',[], @islogical);
            args.addParameter('GroupingVariables', 'Firm', @ischar);
            args.addParameter('weights', [], @ischar);
            args.parse(varargin{:});
            if ~strcmpi(args.Results.GroupingVariables, 'Firm')
                varnames = {args.Results.GroupingVariables, 'Marketid'};
                var = obj.demand.data(:, varnames);
            else
                varnames = {'Firm', 'Marketid'};
                var = table(obj.firm, obj.marketid, 'VariableNames', varnames);
            end
            
            %%%%%% Specific code
            if isempty(args.Results.selection)
                selection = ~isnan(obj2.p);
            else
                selection = args.Results.selection & ~isnan(obj2.p);
            end
            tableCols = {'Costs', 'Price1', 'Price2', 'PriceCh'};
            priceChange = (obj2.p - obj.p ) ./ obj.p;
            res = [var, array2table([ obj.c, obj.p,  obj2.p, priceChange],...
                'VariableNames', tableCols)];            
%                 apc = priceChange' * obj.s / sum(obj.s); 
 %           res.Properties.VariableNames = [varnames, tableCols];
            %%%%%% End specific code
            
            % weights can be set manually
            if ~obj.settings.weightedAverages
                weights = [];
            else
                weights = obj.getMarketShares(args.Results.weights);
            end
            if ~all(selection)
                res = res(selection, :);
                if ~isempty(weights)
                    weights = weights(selection, :);
                end
            end
            x = Estimate.means(res, tableCols, varnames, weights);
            pricetab = Estimate.means(x, tableCols, args.Results.GroupingVariables);
            % Shares
%             varname = 'firm';
%             tab = market.T(:,{varname});
%             tab = [tab, array2table([ market.summary().MarketSh,  market2.summary().MarketSh, ...
%                 market2.summary().MarketSh - market.summary().MarketSh])];
%             tab =  varfun(@sum, tab, 'GroupingVariables', varname);
%             tab(:, {varname, 'GroupCount'}) = [];
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
        
		function f = focNum(obj,  P)
			S = obj.demand.shares( P );
            func = @(p)(obj.demand.shares(p));
            jac = jacobian(func, P);
			f = ( obj.RR .* jac ) * (P - obj.c) + S;
		end 

		function f = margins(obj,  P)
			S = obj.demand.shares( P );
			f = linsolve( (-obj.RR) .* (obj.demand.shareJacobian( S , P)) , S );
        end
        
        function f = fixedPoint(obj, maxit)
			convergence = 1 ;
			P = obj.p0;
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
			obj.s = S;
			f = [i,convergence];
        end
        
        function newmarket = clone(obj)
            newmarket = copy(obj);
        end      
    end
end

