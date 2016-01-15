classdef SimMarket < matlab.mixin.Copyable
    % SimMarket creates simulated data and estimates on it
    % MONTE-CARLO NESTED AND MIXED LOGIT MODELS
    % Simulate data based on model and estimate
    % Products exist with exog probability and count instruments based on the
    % number of products are used to handle endogeneity
    
    properties
        model
        data
        demand
        market
    end
    properties (SetAccess = protected, Hidden = true )
        simDemand % Class with settings for data generation
        estDemand % Settings for estimation
        epsilon
    end
    
    methods
        function obj = SimMarket(varargin)
            % Create model and demand. Alternatively model can be provided:
            % generate model with empty constructor, modify and provide as
            % argument to constructor.
            obj.randdraws();
            if nargin > 0
                obj.model = varargin{1};
            else
                defmodel.endog = false;      % Endog with count instruments or no endog
                defmodel.randomProducts = false; % Let the number of products be random
                defmodel.pricesFromCosts = true; % Simulate price
                defmodel.markets = 100;
                defmodel.products = 5;
                defmodel.types = [];
                defmodel.firm = [];
                
                % p and x mean and sd:
                defmodel.beta = [ 1, 0];
                defmodel.x = [5, 0];
                defmodel.x_vcv = [1, 1];
                defmodel.c = 4;
                defmodel.c_vcv = 1;
                defmodel.gamma = 0;
                
                % Individual and product level shocks:
                defmodel.epsilon_sigma = .1;
                defmodel.sigma_xi = .1;
                
                defmodel.endog_sigma = 0.1; % Degree of corr between x and epsilon
                defmodel.prob_prod = .8;    % Prob of product existing in market
                obj.model = SettingsClass(defmodel);
            end
        end
        
        function rt = create(obj)
            obj.init();
            if obj.model.pricesFromCosts && obj.model.endog
                rt = obj.simulateDemand();
            else
                rt = obj.calculateDemand();
            end
        end
        
        function mr = findCosts(obj, varargin)
            if nargin == 2 
                obj.market.demand = varargin{1};
            else
                obj.market.demand = obj.estDemand;
            end
            if isempty(obj.model.firm)
                obj.market.var.firm = 'productid'; % One product firms
            else
                obj.market.var.firm = 'firm'; % Firm variable has been created
            end
            obj.market.findCosts();
            mr = obj.means( {'q','c',  'p'}, 'productid') ;
            cr = rowfun(@(x,y)(y-x)/y, mr(:,2:3),'OutputVariableNames','Markup');
            mr = [mr, cr];
        end
        
        function results = estimate(obj)
            result = obj.estDemand.estimate();
            if isa(obj.demand, 'RCDemand')
                beta = [-obj.demand.alpha; obj.model.beta'];
            else
                beta = [-obj.demand.alpha; reshape(obj.demand.sigma, [], 1);...
                    obj.model.beta'];
            end
            if strcmpi(obj.estDemand.settings.paneltype, 'fe')
                beta = beta(1:end-1);
            end
            if isa(obj.demand, 'RCDemand')
                truevals = table([beta; obj.demand.sigma]);
            else
                truevals = table(beta);
            end
            truevals.Properties.VariableNames = {'True_val'};
            results = [truevals, result];
        end
        
        function R = means(obj, cols, index)
            R = array2table(splitapply(@mean, obj.data{:, cols}, ...
                obj.data.(index)), 'VariableNames', cols);
        end

    end
    methods (Access = protected, Hidden = true )
        function init(obj)
            % Simulation of dataset
            obj.demand.var.market = 'marketid';
            obj.demand.var.panel = 'productid';
            obj.demand.var.price = 'p';
            obj.demand.var.exog = 'x';
            obj.demand.var.quantity = 'q';
            obj.demand.var.marketsize = 'constant';
            
            obj.simDemand = copy(obj.demand);
            obj.estDemand = copy(obj.demand);
            if isa(obj.demand, 'RCDemand')
                obj.estDemand.sigma = [];
            end
            obj.model.beta = reshape(obj.model.beta, 1, length(obj.model.beta));
            
            obj.createData();
            obj.simDemand.data = obj.data;
            obj.simDemand.d = obj.data.d; % <<<<<<<<<<< init is required
            if isempty(obj.market)
                obj.market = Market();
            end
        end
        
        function createData(obj)
            % createData creates firm, type, and demand and cost shifters.
            % Price is created i
            % Note that demand class is not involved in the creation
            function nt  = nestTree( types )
                nt = (1:types(end))';
                for i = (length(types)-1):-1:1
                    nt = [ reshape(repmat(1:types(i), size(nt, 1), 1), [], 1), ...
                        repmat(nt, types(i), 1)];
                end
            end
            function typemat = createTypes(types, products, markets)
                ntree = nestTree(types);
                typelist = repmat(ntree, ceil(products / size(ntree,2)), 1);
                type = repmat(typelist(1:products, :), markets, 1);
                typemat = array2table(type);
            end
            obj.data = table();
            n = obj.model.markets * obj.model.products;
            obj.data.marketid = reshape(repmat(1:obj.model.markets, ...
                obj.model.products, 1), n, 1);
            obj.data.productid = repmat((1:obj.model.products)', ...
                obj.model.markets, 1);
            if ~isempty(obj.model.types)
                obj.data = [obj.data, createTypes(obj.model.types, ...
                    obj.model.products, obj.model.markets)];
            end
            if ~isempty(obj.model.firm)
                obj.model.firm = reshape(obj.model.firm, length(obj.model.firm),1);
                if length(obj.model.firm) ~= obj.model.products
                    error('The list of firms does not match the number of products')
                end
                obj.data.firm = repmat(obj.model.firm, obj.model.markets, 1);
            end
            
            % epsilon_jt = varepsilon_jt + xi_j
            obj.epsilon =  randn(n, 1) * obj.model.epsilon_sigma + ...
                repmat(randn(obj.model.products, 1) * obj.model.sigma_xi, ...
                obj.model.markets, 1);
            
            % Create random x var:
%            obj.data.x = obj.model.x(2) + randn(n, 1) * obj.model.x_vcv(2);
            x = mvnrnd(obj.model.x, obj.model.x_vcv, n);
            % Create price var to satisfy demand.initAdditional
            obj.data.p = x(:, 1);
            obj.data.x = x(:, 2:end);
            obj.data.constant = ones(n, 1);
            % This should be done in the demand classes:            <<<<<<<<<<<<<<<<
            x0 = [table2array(obj.data(:, 'x')), obj.data.constant];
            obj.data.d = x0 * obj.model.beta' + obj.epsilon;
            
            if obj.model.gamma ~= 0 % Otherwise existing tests fail
                obj.data.w = randn(n, 1);
                obj.data.c = obj.data.constant * obj.model.c ...
                    + obj.model.gamma * obj.data.w ...
                    + randn(n, 1) * obj.model.c_vcv;
            else
                obj.data.c = obj.data.constant * obj.model.c ...
                    + randn(n, 1) * obj.model.c_vcv;
            end
            
            % Random selection of products
            if obj.model.randomProducts
                % probability prob_prod of a product existing in a period
                prodsel = logical(binornd(1, obj.model.prob_prod, n, 1));
                obj.data = obj.data(prodsel, :);
                % TODO: Create count instruments for nests
                instruments = Estimate.countInstruments(obj.data, 'marketid', []);
                obj.data.nprod = instruments(:, 1);
                obj.data.nprod2 = obj.data.nprod .^ 2;
            end
        end
        
        % Function can reset data and i
        function mr = calculateDemand(obj)
            % Create random price:
            n = size(obj.data, 1);
            if obj.model.endog
                % Instruments a la Nevo:
                ninstr = 6;
                A = 2 + 0.2 * randn(n, ninstr);
                M = eye(ninstr)*0.2 + ones(ninstr, ninstr)*0.8;
                inst = A * chol(M);
                % Keep expected effect of instruments on p to be zero
                suminst = sum(inst, 2) - mean(sum(inst, 2));
                obj.data.p = obj.data.p ...
                    + suminst + obj.model.endog_sigma * obj.epsilon(1:n);
            else
                inst = [];
            end
            obj.data = [obj.data, array2table(inst)];
            
            obj.simDemand.init();
            obj.data.q = obj.simDemand.getDemand(obj.data.p);
            
            mr = obj.means({'p', 'q'}, 'productid') ;
            
            display 'Average sum shares'
            disp(mean(accumarray(obj.data.marketid, obj.data.q)))
            if obj.model.endog
                obj.estDemand.var.instruments = sprintf('inst%d ', 1:6);
            end
            obj.estDemand.data = obj.data;
            obj.simDemand.data = obj.data;
        end
        
        function mr = simulateDemand(obj, varargin)
            obj.simDemand.init();
            if nargin > 1
                obj.market.demand = varargin{1};
            else
                obj.market.demand = obj.simDemand;
            end
            if isempty(obj.model.firm)
                obj.market.var.firm = 'productid'; % One product firms
            else
                obj.market.var.firm = 'firm'; % Firm variable has been created
            end
            % Existence of var.costs name could be alternative to:
            obj.market.c = obj.data.c;
            
            obj.market.equilibrium();
            
            obj.data.p = obj.market.p;
            obj.data.q = obj.simDemand.getDemand(obj.data.p);
            
            mr = obj.means({'p', 'q'}, 'productid') ;
            if obj.model.endog && isempty(obj.estDemand.var.instruments)
                if obj.model.randomProducts
                    obj.estDemand.var.instruments = 'nprod nprod2';
                elseif obj.model.gamma ~= 0
                    obj.estDemand.var.instruments = 'w';
                else
                    obj.estDemand.var.instruments = 'c';
                end
            end
            obj.estDemand.data = obj.data;
        end
        
        function randdraws(obj)
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed', 9));
        end
        
        function cpObj = copyElement(obj)
        % Overrides copyElement method:
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the DeepCp object
            cpObj.model = copy(obj.model);
        end

    end
    methods(Static)
        function vals = testEqual(x,y,z)
            diff = abs((y - x)/y);
            vals = [y, x, diff];
            if diff > z
                disp(vals)
                error('True, calculated and percentage diff are greater than %f', z);
            end
        end
    end
end

