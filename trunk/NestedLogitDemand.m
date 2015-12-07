classdef NestedLogitDemand < Estimate
    % Estimation and simulation code for nested logit.
    %   $Id: NestedLogitDemand.m 141 2015-10-08 14:10:51Z d3687-mb $
    
    properties
        alpha
        sigma
        nestlist
        share
        p % Used for CES and simulation
        q 
        ms
        sim % Simulation data
        dummarket 
        d % Utility without the price effect 
        xi % Unobservable utility, residual
    end
    properties (SetAccess = private )
        G % Group membership, each column a group with 1 if member
        H % Subgroup membership
        GG % Group membership 
        HH % Subgroup membership 
        GH % Binary  of subgroup membership in group
    end
    
    methods
        function q = getDemand(obj, p)
            q = zeros(size(obj.dummarket, 1), 1);
            for t = 1:size(obj.dummarket, 2)
                obj.initSimulation(t);
                q(obj.dummarket(:, t)) = obj.shares(p(obj.dummarket(:, t)));
            end            
            if obj.settings.ces
                q = q ./ p;
            end
        end   
        
        function initSimulation(obj, market)
            % initSimulation does the per market initialization
            % Invoke init() for general initialization
            % initSimulation sets obj.d and nesting vectors and matrices
            obj.sim.selection = obj.dummarket(:, market);
            obj.sim.d  = obj.d(obj.dummarket(:, market) );
            obj.sim.market = market;
            if ~isempty(obj.nestlist)
                [~,~, nest1] = unique(obj.data{: , obj.nestlist(1)});
                [~,~, nest2] = unique(obj.data(:, strsplit(strtrim(obj.var.nests))), 'rows');
                if nargin > 0
                    nest1 = nest1(obj.dummarket(:, market));
                    nest2 = nest2(obj.dummarket(:, market));
                end
                obj.G = dummyvar(grp2idx(nest1))';
                obj.GG = obj.G' * obj.G;
                if length(obj.nestlist) == 2
                    obj.H = dummyvar(nest2)';
                    obj.HH = obj.H'*obj.H;
                    obj.GH =( (obj.G*obj.H') > 0);
                    obj.GH = obj.G*obj.H';
                    obj.GH(obj.GH >0)=1;
                end
            end
        end

        function R = estimate(obj, varargin)
            R = estimate@Estimate(obj, varargin{:});
            if ~isempty(obj.results) && isfield(obj.results,'estimate') && ...
                    ~isempty(obj.results.estimate)
                price = obj.X(:, 1);
                priceName = obj.getPriceName();
                if isempty(obj.nestlist)
                    est = obj.results.estimate{priceName, 1}';
                    obj.d = obj.share.ls - price * est';
                elseif length(obj.nestlist) == 1
                    est = obj.results.estimate{[priceName, {'lsjg'}], 1}';
                    obj.d = obj.share.ls - [price, obj.share.lsjg] * est';
                elseif length(obj.nestlist) == 2
                    est = obj.results.estimate{[priceName, {'lsjh', 'lshg'}], 1}';
                    obj.d = obj.share.ls - [price, obj.share.lsjh, obj.share.lshg] * est';
                end
                obj.alpha = -est(1);
                if length(est) > 1
                    obj.sigma = est(2:end);
                end
            end
        end

		function s = shares(obj, P, varargin)
            % shares(p, [market_number])
            if obj.settings.ces
                P = log(P);
            end
            delta = obj.sim.d - obj.alpha .* P; % d is delta without price effect
			if length(obj.nestlist) == 2 
                sigma1 = obj.sigma(1);
                sigma2 = obj.sigma(2);
				ev=exp(delta ./ (1-sigma1));
				ighs = (obj.H*ev);
				igh = obj.H' * ighs;
				evv = ighs .^ ((1-sigma1)/(1-sigma2));
				igs = (obj.GH * evv);
				ig = obj.G' * igs;
				it = sum(igs .^(1-sigma2));
				s = ev .* (igh .^ ((sigma2-sigma1)/(1-sigma2)) ) .* ...
                    (ig .^ (-sigma2)) / (1+it);
			elseif length(obj.nestlist) == 1 
                sigma1 = obj.sigma(1);
				ev=exp(delta ./ (1-sigma1));
				igs = (obj.G * ev);
				ig = obj.G' * igs;
				it = sum(igs .^(1-sigma1));
				s = ev .* (ig .^ (-sigma1)) / (1+it);
			else
				ev=exp(delta );
				it = sum(ev);
				s = ev / (1+it);
			end 
		end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
		% Used in Market.findcosts()
        function S = actualDemand(obj)
            S = obj.share.s;	
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
                
        function sj = shareJacobian(obj, P)
            % shareJacobian([P], [market_number])
            % When P is empty shareJacobian is on actual prices and shares
            % Restrict shares to market t:
            if isempty(P)
                S = obj.share.s(obj.sim.selection);
                P = obj.p(obj.sim.selection);
            else
                S = obj.shares(P);
            end
            other_effect =  - S*S';
			if isempty(obj.nestlist) 	
				own_effect = diag(S);
				sj = -obj.alpha *( own_effect + other_effect );  
            else
                sigma1 = obj.sigma(1);
            
				Sg = obj.GG * S;
				own_effect = 1/(1 - sigma1)*diag(S);
				if length(obj.nestlist) == 2 
                    sigma2 = obj.sigma(2);
					gr_effect = -sigma2/(1 - sigma2)*obj.GG .* ((S ./ Sg)*S');
					Sgh = obj.HH * S;
					gr_effect = gr_effect - (1/(1-sigma1) - 1/(1-sigma2))* ...
                        obj.HH .* ((S ./ Sgh)*S');
				else
					gr_effect = -sigma1/(1 - sigma1)*obj.GG .* ((S ./ Sg)*S');
				end 
				sj = -obj.alpha*( own_effect + gr_effect + other_effect );  
			end
            if obj.settings.ces
                sj = (sj - diag(S) ) * diag(1./P) ;
            end
        end
        
        function [elas, varargout] = elasticities(obj, P)
            function s = sumstats( j, E)
                e = j .* E;
                e(e==0)=[];
                s = [mean(e) std(e) min(e) max(e)];
            end
            s = obj.shares(P);
            n = length(s);
            D = obj.shareJacobian(P)';
            E = diag(P) * D * diag( 1 ./ s );
            elas = [sumstats(eye(n), E)];
            if length(obj.nestlist) == 0
                elas = [elas; sumstats(1 - eye(n), E)]
                rowtit = {'e_ii', 'e_ij'};
            elseif length(obj.nestlist) == 1
                elas = [elas; ...
                    sumstats(obj.GG - eye(n), E); ...
                    sumstats(1 - obj.GG, E)];
                rowtit = {'e_ii', 'e_ij', 'e_ik'};
            elseif length(obj.nestlist) == 2
                elas = [elas; ...
                    sumstats(obj.HH - eye(n), E); ...
                    sumstats(obj.GG - obj.HH , E); ...
                    sumstats(1 - obj.GG, E)];
                rowtit = {'e_ii', 'e_ij', 'e_ik', 'e_il'};
            end
            elas = array2table(elas);
            elas.Properties.VariableNames = {'Mean', 'Std', 'Min', 'Max'};
            elas.Properties.RowNames = rowtit;
            if nargout > 0
                varargout{1} = E;
            end
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
        
        function obj = NestedLogitDemand(varargin)
            obj = obj@Estimate(varargin{:});
            obj.var.setParameters({'quantity','price','nests','marketsize'});
            obj.settings.setParameters({'ces'});
            
            obj.settings.estimateMethod = 'gls';
            obj.settings.ces = false;
            
            obj.results.estimateDescription = 'Nested Logit Demand';             
        end
       
      function name = getPriceName(obj)
            if obj.settings.ces
                name = 'lP';
            else
                name = obj.var.price;
            end
        end
                
        % Invoked from init@Estimate to create some logit specific params
        function names = initAdditional(obj, names, selection)
            if ~isempty(obj.var.nests)
                obj.nestlist = strsplit(strtrim(obj.var.nests));
            end
            lsnames = {[],{'lsjg'},{'lsjh', 'lshg'}};
            names.endog = [obj.getPriceName(), lsnames{length(obj.nestlist)+1}, names.endog];
            
            if isempty(obj.var.market) || isempty(obj.var.price) 
                error('Demand.var.market and price must be specified in model');
            end
            if isempty(selection)
                T = obj.data;
            else
                T = obj.data(selection, :);
            end
            [~,~,id] = unique(T{:, strsplit(strtrim(obj.var.market))}, 'rows');
            obj.marketid = id;
            obj.dummarket = logical(dummyvar(obj.marketid));
            obj.p = T{:, obj.var.price}; % Used in simulation 
            
            % quantity is empty for simulated market
            if ~isempty(obj.var.quantity) && obj.isvar(obj.var.quantity, obj.data)
                obj.var.depvar = 'ls';
                if isempty(obj.var.marketsize) || isempty(obj.var.exog) 
                    error(['Demand.var.exog and marketsize', ...
                        ' must be specified in model']);
                end
                obj.ms = T{: , obj.var.marketsize};
                obj.q = T{: , obj.var.quantity};
                 % init can be invoked several times, so only create once,
                 % unless selection has been reset
                if ~obj.isvar(obj.var.depvar, obj.data) % isempty(obj.share) 
                    obj.share = obj.generateShares(T); 
                    % Not clean !!
                    % obj.data needs shares because depvar uses variable
                    obj.data = [obj.data, obj.share];
                end
                if ~isempty(selection)
                    obj.share = obj.generateShares(T); 
                end
            end
            % This is not clean. Should be added to (renamed) shares struct
            if obj.settings.ces
                obj.data.lP = log(obj.data{:,obj.var.price});
            end
        end
       
        function resultTables(obj)
            resultTables@Estimate(obj);            
            if ~obj.config.quietly
                disp(['Estimate of: ', obj.results.estimateDescription])
                disp(obj.results.estimate);
                s0 = 1 - obj.dummarket' * obj.share.s;
                disp('Share of outside good');
                fprintf('  Mean: %0.3f Min: %0.3f Max: %0.3f \n', ...
                    mean(s0), min(s0), max(s0));
            end            
        end        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods % (Access = protected)
        function S = generateShares(obj, T )
            if obj.settings.ces
                Q = obj.q .* obj.p;
            else
                Q = obj.q;
            end
            subtotal = accumarray(obj.marketid, Q);
            subtotal =  subtotal(obj.marketid,:);
            S = table();
            S.s = Q ./ obj.ms;
            S.s0 = 1 - ( subtotal ./ obj.ms );
            S.ls = log(S.s ./ S.s0);
            market = strsplit(strtrim(obj.var.market));
            if length(obj.nestlist) >= 1
                groupsubtotal =  obj.subtotals(T, Q, ...
                    [market obj.nestlist(1)]);
                S.lsjg = log( Q ./ groupsubtotal);
            end
            if length(obj.nestlist) == 2
                subgroupsubtotal =  obj.subtotals(T, Q, ...
                    [market obj.nestlist]);
                S.lsjh = log( Q ./ subgroupsubtotal);
                S.lshg = log( subgroupsubtotal ./ groupsubtotal );
            end
        end
        
        function f = subtotals(obj, T, sumvar, index )
            %SUBTOTALS Sums sumvar in table T by index category variable list
            %   Mimics Stata by index: egen subtotal = total(sumvar)
                [~, ~, rowIdx] = unique(T( : , index), 'rows');
                subtotal = accumarray(rowIdx, sumvar);
                f =  subtotal(rowIdx,:);
        end
        
        % Standard error of the regression
        % Cameron & Trivedi p 287
        function sd = sdreg(obj)
            sm = obj.share.s - mean(obj.share.s);
            sp = obj.shares(obj.p, 1); % Compute predicted shares
            e = obj.share.s - sp;
            sd = sqrt(e'*e ./(size(obj.X,1)-length(obj.vars2)-length(obj.vars2)));
        end         
    end
end

