classdef Market < Estimate % matlab.mixin.Copyable
    % MARKET Calculates market equilibrium
    %   Based on ownership structure and conduct.
    % $Id: Market.m 134 2015-09-28 19:32:21Z d3687-mb $
    
    properties
        %Varnames
        firm
        q %quantities
        p %prices
        p0 %Initial prices
        c %Costs
        D % Demand class
        
        dampen = 1
        maxit = 3000
        costfunction = 'loglinear'
        estimation = 'ols'
        selection
    end
    properties % (SetAccess = protected )
        s %shares
        R %ownership vector
        RR %ownership tranformation with conduct
        c1, p1, q1 %Views for bootstrapping
        sel % market class with selection
        sim % per market data, set in initSimulation
    end
    
    methods
       function init(obj, varargin)
            if nargin > 1
                selection = varargin{1};
            else
                selection = [];
            end
            obj.D.init();
            % Set unless new firm has been set for example as post merger
            % ownership.
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
%            obj.sim = table(p0, p, q, s, marketid);
            
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
            obj.R = dummyvar(rs)';
            obj.RR  = obj.R' * obj.R;
            if obj.settings.conduct > 0
                obj.RR(obj.RR == 0) = obj.settings.conduct;
            end
        end
        
        function exitflag = equilibrium(obj, varargin)
            obj.init();
            if nargin > 1
                selection = varargin{1};
                obj.selection = selection;
            else
                selection = logical(ones(size(obj.marketid)));
            end
            options = optimoptions(@fsolve,'MaxFunEvals',obj.maxit,'Display', 'off');
            obj.p = zeros(size(obj.marketid));
            obj.s = zeros(size(obj.marketid));
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
                obj.selection = selection;
            else
                selection = logical(ones(size(obj.marketid)));
            end
            obj.c = zeros(size(obj.marketid));
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
			f = ( obj.RR .* obj.D.shareJacobian( P ))*(P - ct)+ S;
		end 

		function f = focNum(obj,  P)
			S = obj.D.shares( P );
            func = @(p)(obj.D.shares(p));
            jac = jacobian(func, P);
			f = ( obj.RR .* jac )*(P - obj.c)+ S;
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
            tableCols = {'Firm', 'Product', 'Costs'};
            demres = obj.D.quantity(obj.p);
            res = table(obj.firm, demres.Product, obj.c, 'VariableNames',tableCols);
            demres.Product = [];
            res = [res, demres];
            if nargin > 1
                sortcols = obj.data(:, varargin{1});
                demres.Price = demres.Price .* demres.MarketSh ;
                demres.Costs = obj.c .* demres.MarketSh ;
                demres = [sortcols, demres];
                res = varfun(@sum, demres, 'GroupingVariables', varargin{1});
                res.GroupCount = [];
                res.Properties.VariableNames = demres.Properties.VariableNames;
                res(:,varargin{1}) = [];
                res.Price = res.Price ./ res.MarketSh ;
                res.Costs = res.Costs ./ res.MarketSh ;
                res.Lerner = (res.Price - res.Costs) ./ res.Price ;
                res(:,{'Quantity', 'Share'}) = [];
            end
        end

        % Collapse mean of table by named column. If variable is RowNames,
        % before invoking create column using:
        % mytab.myvar = mytab.Properties.RowNames.
        function tab = collapseBy(tab, var)
            tab = varfun(@mean, tab, 'GroupingVariables', var);
            tab(:,'GroupCount') = [];
            tab.Properties.VariableNames = tableCols;
            tab(:, var) = [];
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
                
        function [pricetab, apc, varargout] = compare(obj, pn, varargin)
            if nargin == 3 
                varname = varargin(1);
                var = obj.data{:, varname};
            else
                varname = {'Firm'};
                var = obj.firm;
            end
            tableCols = [varname, {'Costs', 'Price1', 'Price2', 'PriceCh'}];
            if isempty(obj.selection)
                priceChange = (pn - obj.p ) ./ obj.p;
                res = table(var, obj.c, obj.p,  pn, ...
                    priceChange, 'VariableNames',tableCols);
                apc = priceChange' * obj.s / sum(obj.s); 
            else
                pn = pn(obj.selection);
                po = obj.p(obj.selection);
                co = obj.c(obj.selection);
                var = var(obj.selection);
                priceChange = (pn - po ) ./ po;
                res = table(var, co, po,  pn, ...
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
        
        function obj = Market(varargin)
            obj = obj@Estimate();
            if nargin > 0
                obj.D = copy(varargin{1});
            end
            varsEstimate = {'market','panel','depvar','exog', ...
                'endog','instruments'};
            obj.var = SettingsClass([varsEstimate, {'firm'}]);
            obj.settings.addprop('conduct');
            obj.settings.conduct = 0;     
        end      
        
        function newmarket = clone(obj)
            newmarket = copy(obj);
        end      
    end
end

