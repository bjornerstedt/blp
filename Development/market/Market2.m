classdef Market2 < Estimate % matlab.mixin.Copyable
    % MARKET Calculates market equilibrium
    %   Based on ownership structure and conduct.
    %   Has separate DemandMarket classes
    % $Id:  $
    
    properties
        %Varnames
        firm
        q %quantities
        p %prices
        p0 %Initial prices
        c %Costs
        % marketid
        D % Demand class
        dm % DemandMarket class
        
        dampen = 1
        maxit = 3000
        costfunction = 'loglinear'
        estimation = 'ols'
    end
    properties % (SetAccess = protected )
        s %shares
        R %ownership vector
        RR %ownership tranformation with conduct
        c1, p1, q1 %Views for bootstrapping
        sel % market class with selection
        isselected = []
        selectindex = []
        selection
    end
    
    methods
        function init(obj)
            obj.period = [];
            if isempty(obj.firm) % Set unless new firm has been set
                obj.firm = obj.D.T{:, obj.var.firm};
            end
            if ~isempty(obj.var.exog)
                init@Estimate(obj);
            end
                obj.initMarkets();
        end
        
        function initSimulation(obj)
            rs = grp2idx(obj.firm);
            obj.R = dummyvar(rs)';
            obj.RR  = obj.R' * obj.R;
            if obj.settings.conduct > 0
                obj.RR(obj.RR == 0) = obj.settings.conduct;
            end
        end
        
        function initMarkets(obj)
            if isempty(obj.period)
                obj.isselected = obj.D.isselected;
                obj.selectindex = find(obj.isselected);
                obj.initPeriods();
            end
        end
        
        function initPeriods(obj)
            periodCopy = copy(obj); 
            obj.period = cell(length(obj.selectindex),1);
            periodCopy.T = [];
            for per = 1:length(obj.selectindex)
                sel = obj.marketid == obj.selectindex(per);
                newperiod = copy(periodCopy); 
                newperiod.T = obj.T(sel, :);
                % Rather ugly code to handle both data creation AND sim.
                if ~isempty(obj.s) 
                    newperiod.s = obj.s(sel,:);
                end
                if ~isempty(obj.p)
                    newperiod.p = obj.p(sel,:);
                end
                if ~isempty(obj.p0)
                    newperiod.p0 = obj.p0(sel,:);
                end
                if ~isempty(obj.c)
                    newperiod.c = obj.c(sel,:);
                end
                newperiod.marketid = obj.marketid(sel,:);
                newperiod.D = obj.D.period{ obj.selectindex(per) };
                newperiod.firm = obj.firm(sel,:);
                newperiod.initSimulation();
                obj.period{per} = newperiod;
            end            
        end
        
        function costs = findCosts(obj)
            if isempty(obj.period)
                % obj.c is used in simulation
                d = obj.D;
                obj.c = obj.p - linsolve(obj.RR .* d.shareJacobian([]), -obj.s );
                costs = obj.c; 
            else
                obj.c = NaN(size(obj.marketid)); 
                for per = 1:length(obj.selectindex)
                    sel = obj.marketid == obj.selectindex(per);
                    costs = findCosts(obj.period{per});
                    obj.c(sel,:) = costs;
                end
            end
        end
        
		function f = equilibrium(obj, varargin)
            if isempty(obj.period)
                options = optimoptions(@fsolve,'MaxFunEvals',obj.maxit,...
                    varargin{:});
                try
                [obj.p,fval,exitflag]  = fsolve(@(x)obj.foc(x), obj.p0, options);
                f = [exitflag];
                catch err
                    warning('The merger simulation did not converge')
                    f = 0;
                end
                obj.results.equilibrium = (f > 0);     
            else
                display('Finding Equilibrium');
                fprintf('Market\tRet. Code\n');
                f = [];
%                 obj.initMarkets(); % Only needed if findCosts not performed
                for per = 1:length(obj.selectindex)
                    index = obj.marketid == obj.selectindex(per);
                    ret = equilibrium(obj.period{per}, 'Display', 'off');
                    obj.p(index,:) = obj.period{per}.p;
                    fprintf('%6d\t%9d\n', per, ret);
                end
            end
            obj.s = obj.D.shares(obj.p);
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
    
		function f = foc(obj,  P)
			S = obj.D.shares( P );
			f = ( obj.RR .* obj.D.shareJacobian( P))*(P - obj.c)+ S;
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
                sortcols = obj.T(:, varargin{1});
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
                var = obj.T{:, varname};
            else
                varname = {'Firm'};
                var = obj.firm;
            end
            tableCols = [varname, {'Costs', 'Price1', 'Price2', 'PriceCh'}];
            priceChange = (pn - obj.p ) ./ obj.p;
            res = table(var, obj.c, obj.p,  pn, ...
                priceChange, 'VariableNames',tableCols);
            if ~isempty(obj.selection)
                res(~obj.selection, :) = [];
            end
            pricetab = varfun(@mean, res, 'GroupingVariables', varname);
            pricetab(:,'GroupCount') = [];
            pricetab.Properties.VariableNames = tableCols;
            pricetab(:,varname) = [];
            apc = priceChange' * obj.s / sum(obj.s); 
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
        
        function obj = Market2(varargin)
            if nargin > 0
                obj.D = varargin{1};
                obj.p0 = obj.D.p;
                obj.p = obj.D.p;
                obj.q = obj.D.q; % Why this duplication?
                obj.s = obj.D.s;   
                obj.marketid = obj.D.marketid;
                obj.T = obj.D.T;
                obj.selection = obj.D.selection;
            end
            obj.var.firm = '';
            obj.settings.conduct = 0;     
        end      
        
        function newmarket = clone(obj)
            newmarket = copy(obj);
        end      
    end
end

% TODO 
% Remove obj.q and obj.s
