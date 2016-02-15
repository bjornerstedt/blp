classdef NLDemand < Estimate
    % Estimation and simulation code for nested logit.
    
    properties
        alpha
        sigma
        d % Utility without the price effect 
        marketid % Protected?
    end
    properties (SetAccess = protected, Hidden = true )
        nestCount = 0
        nest
        share
        p % Used for CES and simulation
        q 
        ms
        sim % Simulation data
        dummarket 
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
            if isempty(obj.alpha)
                obj.alpha = obj.beta(1);
                pars = obj.nestCount + 1;
                if pars > 1
                    obj.sigma = obj.beta(2:pars);
                end
            end
            obj.sim.selection = obj.dummarket(:, market);
            obj.sim.d  = obj.d(obj.dummarket(:, market) );
            obj.sim.market = market;
            if ~isempty(obj.nest)
                obj.G = dummyvar(obj.nest(obj.dummarket(:, market), 1))';
                obj.GG = obj.G' * obj.G;
                if obj.nestCount == 2
                    obj.H = dummyvar(obj.nest(obj.dummarket(:, market), 2))';
                    obj.HH = obj.H' * obj.H;
                    obj.GH =( (obj.G * obj.H') > 0);
                    obj.GH = obj.G * obj.H';
                    obj.GH(obj.GH > 0) = 1;
                end
            end
        end

        function R = estimate(obj, varargin)
        % estimate executes a linear estimate 
        % It uses the settings in NLDemand.settings and NLDemand.var
        
            % If alpha (sigma) and beta have been set, these are
            % displayed as the true values: TODO...
%             if ~isempty(obj.alpha) && ~isempty(obj.sigma)
%             end
            R = estimate@Estimate(obj, varargin{:});
            % These conditions should always hold!:
            if ~isempty(obj.results) && isfield(obj.results,'estimate') && ...
                    ~isempty(obj.results.estimate)
                if obj.settings.ces
                    price = log(obj.p);
                else
                    price = obj.p;
                end
                priceName = obj.getPriceName();
                if isempty(obj.nest)
                    est = obj.results.estimate{priceName, 1}';
                    obj.d = obj.share.ls - price * est';
                elseif obj.nestCount == 1
                    est = obj.results.estimate{[priceName, {'lsjg'}], 1}';
                    obj.d = obj.share.ls - [price, obj.share.lsjg] * est';
                elseif obj.nestCount == 2
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
			if obj.nestCount == 2 
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
			elseif obj.nestCount == 1 
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
            if isempty(obj.share)
                error('Shares have not been calculated')
            end
            S = obj.share.s;	
        end

        function tab = quantity(obj, P)
		% QUANTITY calculates quantities and market shares
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
            % shareJacobian([P])
            % When P is empty shareJacobian is on actual prices and shares
            % Restrict shares to market t:
            if isempty(P)
                S = obj.share.s(obj.sim.selection);
                P = obj.p(obj.sim.selection);
            else
                S = obj.shares(P);
            end
            other_effect =  - S * S';
			if isempty(obj.nest) 	
				own_effect = diag(S);
				sj = -obj.alpha *( own_effect + other_effect );  
            else
                sigma1 = obj.sigma(1);
            
				Sg = obj.GG * S;
				own_effect = 1/(1 - sigma1)*diag(S);
				if obj.nestCount == 2 
                    sigma2 = obj.sigma(2);
					gr_effect = -sigma2/(1-sigma2) * obj.GG .* ((S ./ Sg) * S');
					Sgh = obj.HH * S;
					gr_effect = gr_effect - (1/(1-sigma1) - 1/(1-sigma2))* ...
                        obj.HH .* ((S ./ Sgh) * S');
				else
					gr_effect = -sigma1/(1-sigma1) * obj.GG .* ((S ./ Sg) * S');
				end 
				sj = -obj.alpha*( own_effect + gr_effect + other_effect );  
			end
            if obj.settings.ces
                sj = (sj - diag(S)) * diag(1./P) ;
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
            if obj.nestCount == 0
                elas = [elas; sumstats(1 - eye(n), E)]
                rowtit = {'e_ii', 'e_ij'};
            elseif obj.nestCount == 1
                elas = [elas; ...
                    sumstats(obj.GG - eye(n), E); ...
                    sumstats(1 - obj.GG, E)];
                rowtit = {'e_ii', 'e_ij', 'e_ik'};
            elseif obj.nestCount == 2
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
        
        function [elas] = groupElasticities(obj, P, group)
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
            elas = A * diag(P) * D * A' * diag( 1 ./(A*s) ) ;
            elas = array2table(elas);
            if ~isempty(names)
                elas.Properties.RowNames = names;
                elas.Properties.VariableNames = names;
            end
        end
        
        function obj = NLDemand(varargin)
            obj = obj@Estimate(varargin{:});
            obj.var.setParameters({'quantity','price','nests','marketsize'});
            obj.settings.setParameters({'ces'});
            
            obj.settings.paneltype = 'lsdv';
            obj.settings.estimateMethod = 'gls';
            obj.settings.ces = false;
            
            obj.results.estimateDescription = 'Nested Logit Demand';             
        end
       
        function [xi, beta] = residuals(obj, theta)
            % residuals calculates the residuals for a GMM estimation
            if length(theta) ~= length(obj.sigma) + 1 
                error('Theta is not of the right length');
            end
            obj.alpha = -theta(1);
            obj.sigma = theta(2:end);
            X0 = obj.X(: , (length(theta) + 1):end);
            beta = (X0' * X0) \ (X0' * (obj.y - obj.X(:, 1:length(theta)) * theta));
            beta = [theta; beta];
            xi = obj.y - obj.X * beta;
        end
           
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Access = protected, Hidden = true )
        function name = getPriceName(obj)
            if obj.settings.ces
                name = 'lP';
            else
                name = obj.var.price;
            end
        end
                    
       function  [Xorig, created] = initAdditional(obj)
       % Invoked from init@Estimate to create some logit specific params
       % Creates marketid, dummarket and p
       % If obj.var.quantity is defined and in the data, it creates ms, q,
       % shares, Xorig and y
       % Returns the names of created variables, used in estimation output
            if isempty(obj.var.market) || isempty(obj.var.price) 
                error('Demand.var.market and price must be specified in model');
            end
            if ~isempty(obj.var.nests)
                nestlist = strsplit(strtrim(obj.var.nests));
                obj.nestCount = length(nestlist);
                obj.nest = [];
                for i = 1:length(nestlist)
                    [~,~, nesti] = unique(obj.data{:, nestlist(1:i)}, 'rows');
                    obj.nest = [obj.nest, nesti];
                end
            end
            lsnames = {[], {'lsjg'}, {'lsjh', 'lshg'}};
            if ~isempty(obj.var.quantity) && obj.isvar(obj.var.quantity, obj.data)
                created = [obj.getPriceName(), lsnames{obj.nestCount+1}];
            else
                if isempty(obj.alpha)
                    error('Either quantities or alpha have to be specified')
                end
                created = obj.getPriceName();
            end
            [~,~,id] = unique(obj.data{:, strsplit(strtrim(obj.var.market))}, 'rows');
            obj.marketid = id;
            obj.dummarket = logical(dummyvar(obj.marketid));
            obj.p = obj.data{:, obj.var.price};  
            
            if obj.settings.ces
                Xorig = [log(obj.p)];
            else
                Xorig = [obj.p];
            end
            % quantity is empty for simulated market
            if ~isempty(obj.var.quantity) && obj.isvar(obj.var.quantity, obj.data)
                if isempty(obj.var.marketsize) || isempty(obj.var.exog) 
                    error(['Demand.var.exog and marketsize', ...
                        ' must be specified in model']);
                end
                obj.ms = obj.data{: , obj.var.marketsize};
                obj.q = obj.data{: , obj.var.quantity};
                obj.share = obj.generateShares(obj.data);
                obj.y = obj.share.ls;
                sh = obj.share{:, lsnames{obj.nestCount+1}};
                Xorig = [Xorig, sh];
            end
        end
        
        function wa = useValueShares(obj)
            wa = obj.settings.ces;
        end
        
        function resultTables(obj)
            resultTables@Estimate(obj);  
            ts = obj.dummarket' * obj.share.s;
            obj.results.totalShares = sprintf('  Mean: %0.3f Min: %0.3f Max: %0.3f \n', ...
                mean(ts), min(ts), max(ts));
            if ~obj.config.quietly
                disp(['Estimate of: ', obj.results.estimateDescription])
                disp(obj.results.estimate);
                s0 = 1 - obj.dummarket' * obj.share.s;
                disp('Share of outside good');
                fprintf('  Mean: %0.3f Min: %0.3f Max: %0.3f \n', ...
                    mean(s0), min(s0), max(s0));
            end
        end
        
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
            if obj.nestCount >= 1
                market = strsplit(strtrim(obj.var.market));
                nestlist = strsplit(strtrim(obj.var.nests));
                groupsubtotal =  obj.subtotals(T, Q, ...
                    [market nestlist(1)]);
                S.lsjg = log( Q ./ groupsubtotal);
                
                if obj.nestCount == 2
                    subgroupsubtotal =  obj.subtotals(T, Q, ...
                        [market nestlist]);
                    S.lsjh = log( Q ./ subgroupsubtotal);
                    S.lshg = log( subgroupsubtotal ./ groupsubtotal );
                end
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

