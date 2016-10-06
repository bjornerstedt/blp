classdef NLDemand < Estimate
    % Estimation and simulation code for nested logit.
    
    properties
        alpha
        sigma
        d % Utility without the price effect 
    end
    properties (SetAccess = protected, Hidden = true )
        marketid % Protected?
        nestCount = 0
        nest
        share
        ms
        simMarket % Simulation market number
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
                q(obj.dummarket(:, t)) = obj.shares(p(obj.dummarket(:, t)), t);
            end            
            if obj.settings.ces
                q = q ./ p;
            end
            q = q .* obj.ms;
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
            obj.simMarket = market;
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

        function truevals = setTrueResults(obj, beta)
            beta = [-obj.alpha; reshape(obj.sigma, [], 1); beta'];
            if strcmpi(obj.settings.paneltype, 'fe')
                beta = beta(1:end-1);
            end
            truevals = table(beta);
            truevals.Properties.VariableNames = {'True_val'};
            obj.results.trueValues = truevals;
        end

        function R = estimate(obj, varargin)
            % estimate executes a linear estimate
            % It uses the settings in NLDemand.settings and NLDemand.var
            
            R = estimate@Estimate(obj, varargin{:});

            obj.calibrate();

            if ~isempty(obj.results) && isfield(obj.results,'trueValues')
                R = [obj.results.trueValues, obj.results.estimate];
            end
        end

		function s = calibrate(obj)
            price = obj.data{:, obj.var.price};
            if obj.settings.ces
                price = log(price);
            end
            est = obj.beta(1:(1 + obj.nestCount))';
            if isempty(obj.nest)
                shdata = price ;
            else
                shdata = [price, obj.share{:,(end-obj.nestCount+1):end}];
                obj.results.sigma = obj.sigma;
            end
            obj.d = obj.share.ls - shdata * est';
            obj.alpha = -obj.beta(1);
            obj.results.alpha = obj.alpha;
            if length(est) > 1
                obj.sigma = est(2:end);
            end
        end

		function s = shares(obj, P, market)
            % shares(p, market_number)
            
            % Only initialize if the market has changed:
            if market ~= obj.simMarket
                obj.initSimulation(market)
                obj.simMarket = market;
            end
            if obj.settings.ces
                P = log(P);
            end
            delta = obj.d(obj.dummarket(:, market) ) - obj.alpha .* P; 
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

		% Calculate quantities and market shares
        function tab = quantity(obj, P, t)
            tableCols = {'Product', 'Price', 'Quantity', 'MarketSh', 'Share'};
            s = obj.shares(P, t);
            if obj.settings.ces
                q = s .* obj.ms ./ P;
            else
                q = s .* obj.ms;
            end
            m = q / sum(q); 
            tab = table(obj.panelid, P, q, m, s);
            tab.Properties.VariableNames = tableCols;
        end
                
        function sj = shareJacobian(obj, P, market)
            % shares initializes market if necessary:
            S = obj.shares(P, market);
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
        
        function [elas, varargout] = elasticities(obj, selection, varargin)
            % [table, mat] = elasticities(selection, [price], finds the average
            % elasticities by nesting. The second result is the complete
            % matrix.
            function s = sumstats( E, j)
            % All the elasticities within the same group level are identical
            % along rows. sumstats uses max to sele
                e = j .* E;
                % Because each row in E has a unique value, the following
                % command will find it, whether it is positive or negative.
                e = max(e,[], 2) + min(e,[],2);
                s = [mean(e) std(e) min(e) max(e)];
            end
            function s = avstats( E, j)
            % All the elasticities within the same group level are identical
            % along rows. sumstats uses max to sele
                e = j .* E;
                % Take row means of nonzero elements:
                e = sum(e, 2) ./ sum(e ~= 0,2);
                s = [mean(e) std(e) min(e) max(e)];
            end
            args = inputParser;
            args.addRequired('selection', @islogical);
            args.addParameter('price', obj.data{:, obj.var.price}, @isnumeric);
            args.addParameter('group', [], @ischar);
            args.parse(selection, varargin{:});            
            P = args.Results.price;
            marketId = unique(obj.marketid(selection));
            if length(marketId) >1
                error('Elasticities in multiple markets not supported')
            end
            P = P(obj.dummarket(:, marketId));
            % shares() initializes period, if necessary:
            s = obj.shares(P, marketId);
            n = length(s);
            D = obj.shareJacobian(P, marketId)';
            E = diag(P) * D * diag( 1 ./ s );
            elas = [sumstats(E, eye(n))];
            if ~isempty(args.Results.group)
                G = dummyvar( obj.data{:, args.Results.group} );
                GG = G * G';
                elas = [elas; ...
                    obj.sumstats(E, GG - eye(n)); ...
                    obj.sumstats(E, 1 - GG)];
                rowtit = {'e_ii', 'e_ij', 'e_ik'};
            elseif obj.nestCount == 0
                elas = [elas; sumstats(E, 1 - eye(n))]
                rowtit = {'e_ii', 'e_ji'};
            elseif obj.nestCount == 1
                elas = [elas; ...
                    sumstats(E, obj.GG - eye(n)); ...
                    sumstats(E, 1 - obj.GG)];
                rowtit = {'e_ii', 'e_ji', 'e_ki'};
            elseif obj.nestCount == 2
                elas = [elas; ...
                    sumstats(E, obj.HH - eye(n)); ...
                    sumstats(E, obj.GG - obj.HH ); ...
                    sumstats(E, 1 - obj.GG)];
                rowtit = {'e_ii', 'e_ji', 'e_ki', 'e_ll'};
            end
            elas = array2table(elas);
            elas.Properties.VariableNames = {'Mean', 'Std', 'Min', 'Max'};
            elas.Properties.RowNames = rowtit;
            if nargout > 0
                varargout{1} = E;
            end
        end
        
        function [elas] = groupElasticities(obj, group, selection, varargin)
            if nargin == 4
                P = varargin{1};
            else
                P = obj.data{:, obj.var.price};
            end
            marketId = unique(obj.marketid(selection));
            if length(marketId) >1
                error('Elasticities in multiple markets not supported')
            end
            P = P(obj.dummarket(:, marketId));
            group = obj.data{obj.dummarket(:, marketId), group};
            if iscategorical(group)
                [names,~,group] = unique(group);
                names = matlab.lang.makeValidName(cellstr(char(names)));
            else
                names = [];
            end
            A = dummyvar(group)';
            s = obj.shares(P, marketId);
            D = obj.shareJacobian(P, marketId)';
            elas = A * diag(P) * D * A' * diag( 1 ./(A*s) ) ;
            elas = array2table(elas);
            if ~isempty(names)
                elas.Properties.RowNames = names;
                elas.Properties.VariableNames = names;
            end
        end
        
        function obj = NLDemand(varargin)
            obj = obj@Estimate(varargin{:});
            obj.var.setParameters({'quantity','price','nests','marketSize'});
            obj.settings.setParameters({'ces'});
            
            obj.settings.paneltype = 'lsdv';
            obj.settings.estimateMethod = '2sls';
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
            if obj.settings.ces
                Xorig = [log(obj.data{:, obj.var.price})];
            else
                Xorig = obj.data{:, obj.var.price};
            end
            if isempty(obj.var.marketSize) || isempty(obj.var.exog)
                error(['Demand.var.exog and marketSize', ...
                    ' must be specified in model']);
            end
            obj.ms = obj.data{: , obj.var.marketSize};
            % quantity is empty for simulated market
            if ~isempty(obj.var.quantity) && obj.isvar(obj.var.quantity, obj.data)
                obj.share = obj.generateShares();
                obj.y = obj.share.ls;
                sh = obj.share{:, lsnames{obj.nestCount+1}};
                Xorig = [Xorig, sh];
                obj.results.s0 = ...
                    mean(accumarray(obj.marketid, obj.share.s0, [], @mean));
            end
            obj.simMarket = 0;
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
        
        function S = generateShares(obj)
            p = obj.data{:, obj.var.price};
            q = obj.data{: , obj.var.quantity};
           if obj.settings.ces
                Q = q .* p;
            else
                Q = q;
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
                groupsubtotal =  obj.subtotals( Q, ...
                    [market nestlist(1)]);
                if obj.nestCount == 1
                    
                    S.lsjg = log( Q ./ groupsubtotal);
                elseif obj.nestCount == 2
                    subgroupsubtotal =  obj.subtotals( Q, ...
                        [market nestlist]);
                    S.lsjh = log( Q ./ subgroupsubtotal);
                    S.lshg = log( subgroupsubtotal ./ groupsubtotal );
                end
            end
        end
        
        function f = subtotals(obj, sumvar, index )
            %SUBTOTALS Sums sumvar in data by index category variable list
            %   Mimics Stata by index: egen subtotal = total(sumvar)
                [~, ~, rowIdx] = unique(obj.data( : , index), 'rows');
                subtotal = accumarray(rowIdx, sumvar);
                f =  subtotal(rowIdx,:);
        end
        
        % Standard error of the regression
        % Cameron & Trivedi p 287
        function sd = sdreg(obj)
        % sdreg is not currently used. Note that to use, the p variable has
        % to be restricted to selection
            p = obj.data{:, obj.var.price};
            sm = obj.share.s - mean(obj.share.s);
            sp = obj.shares(p, 1); % Compute predicted shares
            e = obj.share.s - sp;
            sd = sqrt(e'*e ./(size(obj.X,1)-length(obj.vars2)-length(obj.vars2)));
        end   
        
    end
end

