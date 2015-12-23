classdef Market < Estimate % matlab.mixin.Copyable
    % MARKET Calculates market equilibrium
    %   Based on ownership structure and conduct.
    
    properties
        firm
        q %quantities
        p %prices
        p0 %Initial prices
        c %Costs
        D % Demand class        
    end
    properties (SetAccess = protected, Hidden = true )
        maxit = 3000
        s %shares
        RR %ownership tranformation with conduct
        c1, p1, q1 %Views for bootstrapping
        sel % market class with selection
        sim % per market data, set in initSimulation
    end
    
    methods
       function init(obj, varargin)
            args = inputParser;
            args.addParameter('selection', [],  @islogical );
            args.parse(varargin{:});
            selection = args.Results.selection;
            obj.D.init(); % init selection??
            % Set unless new firm has been set for example as post merger
            % ownership.
            if isempty(obj.var.firm) 
                error('var.firm has to be specified');
            end
            if isempty(obj.firm)
                obj.firm = obj.D.data{:, obj.var.firm};
            end

            % Market 'inherits' variables from Demand, to avoid unnecessary
            % specification
            obj.data = obj.D.data;
            % This is for cost estimation:
            if ~isempty(obj.var.exog)
                init@Estimate(obj);
            end
            % Use weighting of averages as in demand
            obj.settings.weights = obj.D.settings.weights;
            % Can probably join these in a table directly, for sim
            obj.p = obj.D.p;
            obj.q = obj.D.q; 
            % With simulated demand, obj.share has not been set
            if ~isempty(obj.D.share)
                obj.p0 = obj.D.p;
                obj.s = obj.D.share.s;   % alt: obj.D.actualDemand()
            else
                obj.p0 = obj.c;                
            end
            obj.marketid = obj.D.marketid;
            obj.panelid = obj.D.panelid;            
            if ~isempty(selection) 
                if isempty(obj.D.selection)  
                    % Rather ugly code to handle both data creation AND 
                    if ~isempty(obj.q) 
                        obj.q = obj.q(selection,:);
                    end
                    if ~isempty(obj.p)
                        obj.p = obj.p(selection,:);
                    end
                    if ~isempty(obj.p0)
                        obj.p0 = obj.p0(selection,:);
                    end
                    if ~isempty(obj.c)
                        obj.c = obj.c(selection,:);
                    end
                    obj.s = obj.D.share.s(selection,:);
                end
                if length(obj.firm) ~= length(obj.s)
                    obj.firm = obj.firm(selection,:);
                end
            end
        end
        
        function initSimulation(obj, market_number)
            obj.D.initSimulation(market_number);
            rs = obj.firm(obj.marketid == market_number);
            %ownership vector:
            R = dummyvar(rs)';
            obj.RR  = R' * R;
            if obj.settings.conduct > 0
                obj.RR(obj.RR == 0) = obj.settings.conduct;
            end
        end
        
        function exitflag = equilibrium(obj, varargin)
            obj.init();
            if nargin > 1
                selection = varargin{1};
            else
                selection = logical(ones(size(obj.marketid)));
            end
            selection = selection & ~isnan(obj.c);
            options = optimoptions(@fsolve,'MaxFunEvals',obj.maxit,'Display', 'off');
            obj.p = nan(size(obj.marketid));
            obj.s = nan(size(obj.marketid));
            marketid_list = unique(obj.marketid(selection));
            convcount = 0;
            for i = 1:length(marketid_list)
                selection = obj.marketid == marketid_list(i);
                t = marketid_list(i);
                obj.initSimulation(t);
                %                try
                [pt,~,exitflag] = fsolve( ...
                    @(x)obj.foc(x, obj.c(obj.marketid == t)), ...
                    obj.p0(selection), options);
                %                 catch err
                %                     warning('The merger simulation did not converge')
                %                     exitflag = 0;
                %                 end
                obj.p(selection) = pt;
                obj.s(selection) = obj.D.shares(pt);
                obj.results.equilibrium = (exitflag > 0);
                if exitflag == 1
                    convcount = convcount+1;
                end
            end
            display(sprintf('Simulation converged for %d of %d markets', ...
                convcount, length(marketid_list)));
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
				P = obj.dampen*Pn + (1 - obj.dampen)*P;
			end 
			if i >= maxit 
				warning( 'Max number of iterations exceeded.' );
				convergence = 0;
			end 
			S = obj.D.shares( P ); %This function should be called demand
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
    
        function findCosts(obj, varargin)
            obj.init();
            if nargin > 1
                selection = varargin{1};
            else
                selection = logical(ones(size(obj.marketid)));
            end
            obj.c = nan(size(obj.marketid));
            marketid_list = unique(obj.marketid(selection));
            for i = 1:length(marketid_list)
                selection = obj.marketid == marketid_list(i);
                t = marketid_list(i);
                obj.initSimulation(t);
                ct = obj.p(selection) - ...
                    linsolve( obj.RR .* obj.D.shareJacobian([]), ...
                    -obj.s(selection) );
                obj.c(selection) = ct;
            end
        end

        function theta = gmm_estimate(obj, theta)
            % Simultaneous 2SLS estimate of demand and costs over alpha
            
            WC = inv(obj.X' * obj.X);
            W = inv(obj.D.Z' * obj.D.Z);
            options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'MaxIter',50);
            
            [theta] = fminunc(@(x)objective(x, obj), theta, options);
            [~, beta ] = obj.D.residuals(theta);
            [~, gamma ] = obj.residuals();
            theta = [beta; gamma];
            
            function val = objective(theta, obj)
                % Residuals as function of demand params:
                xi = obj.D.residuals(theta);
                obj.findCosts();
                eta = obj.residuals();
                
                xiZ = xi' * obj.D.Z;
                etaZ = eta' * obj.X;
                val = xiZ * W * xiZ' + etaZ * WC * etaZ';
            end
        end
        
        function [xi, gamma] = residuals(obj)
            gamma = (obj.X' * obj.X) \ ( obj.X' * obj.c);
            xi = obj.c - obj.X * gamma;
        end
        
		function f = foc(obj,  P, ct)
			S = obj.D.shares( P );
			f = ( obj.RR .* obj.D.shareJacobian( P )) * (P - ct) + S;
		end 

		function f = focNum(obj,  P)
			S = obj.D.shares( P );
            func = @(p)(obj.D.shares(p));
            jac = jacobian(func, P);
			f = ( obj.RR .* jac ) * (P - obj.c) + S;
		end 

		function f = margins(obj,  P)
			S = obj.D.shares( P );
			f = linsolve( (-obj.RR) .* (obj.D.shareJacobian( S , P)) , S );
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
        
        function res = summary(obj, varargin)
            % SUMMARY(group_vars) creates summary by varargin variable
            % Mean over markets of: p, c, Lerner index and market share.
            % Within group, means can be unweighted or weighted by quantity
            % or by value. Without group_vars, reporting is at the product level.
            % Output is as a table.
            res = array2table([obj.marketid, obj.p, obj.c]);
            % Calculate market shares by quantity or value
            if obj.settings.valueShares % Note that valueShares should be true for ces.
                qu = obj.q * obj.p;
            else
                qu = obj.q;
            end
            qbar = accumarray(obj.panelid, qu, [], @sum);
            qbar = qbar(obj.panelid, :);
            res.m = qu ./ qbar;
            res.Properties.VariableNames = {'Price', 'Costs', 'MarketSh', 'Market'};
            res.Lerner = (res.Price - res.Costs) ./ res.Price ;
            if nargin > 1
                groupvar = varargin{1};
                [~,~, sortcols] = unique(obj.data(:, varargin));
                res.(varargin{1}) = sortcols;
            else
                groupvar = 'Firm';
                res.Firm = obj.firm;
            end
            if obj.settings.weightedAverages
                res.Price = res.Price .* res.MarketSh ;
                res.Costs = res.Costs .* res.MarketSh ;
                res = varfun(@sum, res, 'GroupingVariables', [{groupvar}, 'Market']);
                res.GroupCount = [];
                res.Properties.VariableNames = {groupvar, 'Market', 'Price', 'Costs', 'MarketSh', 'Lerner'};
                res.Price = res.Price ./ res.MarketSh ;
                res.Costs = res.Costs ./ res.MarketSh ;
                res.Market = [];
            end
            res = varfun(@mean, res, 'GroupingVariables', groupvar);
            res.GroupCount = [];
            res.Properties.VariableNames = {groupvar,  'Price', 'Costs', 'MarketSh', 'Lerner'};
            res.MarketSh = [];
        end

        function tab = compareAll(obj, other)
            tableCols = {'Firm', 'Product'};
            obj1 = obj.summary();
            obj2 = other.summary();
            firstCols = obj1(:,tableCols);
            obj1(:,tableCols) = [];
            obj2(:,tableCols) = [];
            varnames = obj1.Properties.VariableNames;
            obj1 = table2array(obj1);
            obj2 = table2array(obj2);
            diff = array2table((obj2 - obj1) ./ obj1);
            diff.Properties.VariableNames = ...
                cellfun(@(x) {sprintf('%sCh', x)}, varnames );
            diff.QuantityCh = []; % Same as market share change
            tab = [firstCols, diff];
        end
                
        function [pricetab, apc, varargout] = compare(obj, obj2, varargin)
            args = inputParser;
  %          args.addRequired('obj2', @(x)isa(x, 'NestedLogitDemand'));
            args.addParameter('selection',[], @islogical);
            args.addParameter('groupvar', 'Firm', @ischar);
            args.parse(varargin{:});
            if ~strcmpi(args.Results.groupvar, 'Firm')
                varname = args.Results.groupvar;
                var = obj.data{:, varname};
            else
                varname = {'Firm'};
                var = obj.firm;
            end
            if isempty(args.Results.selection)
                selection = ~isnan(obj2.p);
            else
                selection = args.Results.selection & ~isnan(obj2.p);
            end
            tableCols = [varname, {'Costs', 'Price1', 'Price2', 'PriceCh'}];
            if all(selection)
                priceChange = (obj2.p - obj.p ) ./ obj.p;
                res = table(var, obj.c, obj.p,  obj2.p, ...
                    priceChange, 'VariableNames', tableCols);
                apc = priceChange' * obj.s / sum(obj.s); 
            else
                p2 = obj2.p(selection);
                po = obj.p(selection);
                co = obj.c(selection);
                var = var(selection);
                priceChange = (p2 - po ) ./ po;
                res = table(var, co, po,  p2, ...
                    priceChange, 'VariableNames',tableCols);            
            end
            pricetab = varfun(@mean, res, 'GroupingVariables', varname);
            pricetab(:,'GroupCount') = [];
            pricetab.Properties.VariableNames = tableCols;
            pricetab(:,varname) = [];
            % Shares
%             varname = 'firm';
%             tab = market.T(:,{varname});
%             tab = [tab, array2table([ market.summary().MarketSh,  market2.summary().MarketSh, ...
%                 market2.summary().MarketSh - market.summary().MarketSh])];
%             tab =  varfun(@sum, tab, 'GroupingVariables', varname);
%             tab(:, {varname, 'GroupCount'}) = [];
            
            if nargout == 3
                varargout{1} = priceChange;
            end
        end
        
        function [pricetab, apc, varargout] = compare2(obj, obj2, varargin)
            args = inputParser;
  %          args.addRequired('obj2', @(x)isa(x, 'NestedLogitDemand'));
            args.addParameter('selection',[], @islogical);
            args.addParameter('GroupingVariables', 'Firm', @ischar);
            args.addParameter('weights', [], @ischar);
            args.parse(varargin{:});
            if ~strcmpi(args.Results.GroupingVariables, 'Firm')
                varnames = {args.Results.GroupingVariables, obj.var.market};
                var = obj.data{:, varnames};
            else
                varnames = {'Firm', obj.var.market};
                var = [obj.firm, obj.marketid];
            end
            if isempty(args.Results.selection)
                selection = ~isnan(obj2.p);
            else
                selection = args.Results.selection & ~isnan(obj2.p);
            end
            tableCols = [varnames, {'Costs', 'Price1', 'Price2', 'PriceCh'}];
            if all(selection)
                priceChange = (obj2.p - obj.p ) ./ obj.p;
                res = table(var, obj.c, obj.p,  obj2.p, ...
                    priceChange, 'VariableNames', tableCols);
                apc = priceChange' * obj.s / sum(obj.s); 
            else
                p2 = obj2.p(selection);
                po = obj.p(selection);
                co = obj.c(selection);
                var = var(selection);
                priceChange = (p2 - po ) ./ po;
                res = table(var, co, po,  p2, ...
                    priceChange, 'VariableNames',tableCols);            
            end
            if isempty(arg.Results.weights)
                weights = obj.settings.weights;
            else
                weights = arg.Results.weight;
            end
            res.(weights) = obj.data.(weights);
            pricetab = obj.summarise(@mean, res, 'GroupingVariables', ...
                varnames, 'InputVariables', tableCols(2:end), ...
                'weights', weights);
            
            % Shares
%             varname = 'firm';
%             tab = market.T(:,{varname});
%             tab = [tab, array2table([ market.summary().MarketSh,  market2.summary().MarketSh, ...
%                 market2.summary().MarketSh - market.summary().MarketSh])];
%             tab =  varfun(@sum, tab, 'GroupingVariables', varname);
%             tab(:, {varname, 'GroupCount'}) = [];
            
            if nargout == 3
                varargout{1} = priceChange;
            end
        end
        
        function obj = Market(varargin)
            obj = obj@Estimate();
            if nargin > 0
                obj.D = copy(varargin{1});
            end
            varsEstimate = {'market','panel','depvar','exog', ...
                'endog','instruments'};
            obj.var = SettingsClass([varsEstimate, {'firm'}]);
            obj.settings.setParameters({'conduct', 'dampen', ...
                'weightedAverages', 'valueShares'});
            obj.settings.conduct = 0;  
            obj.settings.dampen = 1;
            obj.settings.weightedAverages = true;
            obj.settings.valueShares = false;
        end      
        
        function newmarket = clone(obj)
            newmarket = copy(obj);
        end      
    end
end

