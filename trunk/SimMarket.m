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
        epsilon
    end
    
    methods
        function obj = SimMarket(varargin)
            % Create model and demand. Alternatively model can be provided:
            % generate model with empty constructor, modify and provide as
            % argument to constructor.
            SimMarket.randDraws();
            if nargin > 0
                obj.model = varargin{1};
            else
                defmodel.endog = false;      % Endog with count instruments or no endog
                defmodel.randomProducts = false; % Let the number of products be random
                defmodel.pricesFromCosts = true; % Simulate price
                defmodel.markets = 100;
                defmodel.products = 5;
                defmodel.types = [];
                defmodel.typeList = [];
                defmodel.firm = [];
                
% beta does not include alpha, but includes constant:
                defmodel.beta = [ 1, 0];
% x is the expected values, const not included. The value for p only
% matters if market is calculated, not simulated. 
                defmodel.x = [5, 0];
                defmodel.x_vcv = [1, 1];

                defmodel.gamma = 0;
                defmodel.c = 4; % Should be part of gamma
                defmodel.w = 0;
                defmodel.w_vcv = 1;
                defmodel.eta = 1; % Cost func error term
                % No fixed effects
                
                % Individual and product level shocks:
                defmodel.epsilon = .1;
                defmodel.xi = .1;
                defmodel.varepsilon = 0; % Utility variation not observed by firms
                
                defmodel.endog_vcv = 0.1; % Degree of corr between x and epsilon
                defmodel.productProbability = .8;    % Prob of product existing in market
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
                obj.market.demand = obj.demand;
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
               
        function R = means(obj, cols, index)
            R = splitapply(@mean, obj.data{:, cols}, ...
                obj.data.(index));
            if length(cols) == size(obj.data{:, cols},2)
                R = array2table(R, 'VariableNames', cols);
            end
        end

    end
    methods (Access = protected, Hidden = true )
        function init(obj)
            % Simulation of dataset
            obj.demand.var.market = 'marketid';
            obj.demand.var.panel = 'productid';
            varNames = {'price', 'quantity', 'exog', 'marketSize'};
            varValues = {'p', 'q', 'x', 'constant'};
            for i = 1:length(varNames)
                if isempty(obj.demand.var.(varNames{i}))
                    obj.demand.var.(varNames{i}) = varValues{i};
                end
            end
            obj.simDemand = copy(obj.demand);
            obj.demand.setTrueResults( obj.model.beta);
            if isa(obj.demand, 'RCDemand')
                % obj.sigma has to be cleared in order for estimation to select
                % a random starting point. Otherwise obj.sigma is used.
                obj.demand.sigma = [];
                obj.simDemand.settings.sigma0 = [];
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
                typeList = repmat(ntree, ceil(products / size(ntree,2)), 1);
                type = repmat(typeList(1:products, :), markets, 1);
                typemat = array2table(type);
            end
            obj.data = table();
            n = obj.model.markets * obj.model.products;
            obj.data.marketid = reshape(repmat(1:obj.model.markets, ...
                obj.model.products, 1), n, 1);
            obj.data.productid = repmat((1:obj.model.products)', ...
                obj.model.markets, 1);
            typeNames = [];
            if ~isempty(obj.model.types)
                typeVars = createTypes(obj.model.types, ...
                    obj.model.products, obj.model.markets);
                obj.data = [obj.data, typeVars];
                typeNames = typeVars.Properties.VariableNames;
            end
            if ~isempty(obj.model.typeList)
                if min(size(obj.model.typeList)) == 1
                    obj.model.typeList = reshape(obj.model.typeList, length(obj.model.typeList),1);
                end
                if size(obj.model.typeList, 1) ~= obj.model.products
                    error('The number of rows types does not match the number of products')
                end
                obj.data.type = repmat(obj.model.typeList, obj.model.markets, 1);
                typeNames = {'type'};
            end
            if ~isempty(obj.model.firm)
                obj.model.firm = reshape(obj.model.firm, length(obj.model.firm),1);
                if length(obj.model.firm) ~= obj.model.products
                    error('The list of firms does not match the number of products')
                end
                obj.data.firm = repmat(obj.model.firm, obj.model.markets, 1);
            end
            
            % epsilon_jt = varepsilon_jt + xi_j
            obj.epsilon =  randn(n, 1) * obj.model.epsilon + ...
                repmat(randn(obj.model.products, 1) * obj.model.xi, ...
                obj.model.markets, 1);
            
            % Create random x var:
            x = mvnrnd(obj.model.x, obj.model.x_vcv, n);
            % Create price var to satisfy demand.initAdditional
            obj.data.p = x(:, 1);
            x = x(:, 2:end);
            obj.data = [obj.data, array2table(x)];
            obj.data.constant = ones(n, 1);
            % This should be done in the demand classes:            <<<<<<<<<<<<<<<<
            x0 = [x, obj.data.constant];
            obj.data.d = x0 * obj.model.beta' + obj.epsilon;
            
            if obj.model.gamma ~= 0 % Otherwise existing tests fail
                obj.data.w = obj.model.w + randn(n, 1) * obj.model.w_vcv;
                obj.data.c = obj.data.constant * obj.model.c ...
                    + obj.model.gamma * obj.data.w ...
                    + randn(n, 1) * obj.model.eta;
            else
                obj.data.c = obj.data.constant * obj.model.c ...
                    + randn(n, 1) * obj.model.eta;
            end
            
            % Random selection of products
            if obj.model.randomProducts
                % probability prob_prod of a product existing in a period
                prodsel = logical(binornd(1, obj.model.productProbability, n, 1));
                obj.data = obj.data(prodsel, :);
                % Remove markets with only one product
                sp = find(sum(dummyvar(obj.data.marketid)) ==1 );
                for i = 1:length(sp)
                    obj.data(obj.data.marketid == sp(i), :) = [];
                end
                % TODO: Create count instruments for nests
                countVars = [];
                if ~isempty(obj.model.firm)
                    countVars = { 'firm' };
                end
                if ~isempty(typeNames)
                    countVars = [countVars, typeNames];
                end
                % xxxx
                instruments = Estimate.countInstruments(obj.data, 'marketid', countVars);
                obj.data.nprod = instruments;
    %            obj.data.nprod2 = obj.data.nprod .^ 2;
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
                    + suminst + obj.model.endog_vcv * obj.epsilon(1:n);
            else
                inst = [];
            end
            obj.data = [obj.data, array2table(inst)];
            if obj.model.varepsilon > 0
                obj.simDemand.d = obj.simDemand.d + ...
                    randn(size(obj.data,1), 1) * obj.model.varepsilon;
            end            
            obj.simDemand.init();
            obj.data.q = obj.simDemand.getDemand(obj.data.p);
            if obj.model.varepsilon > 0
                obj.data.q = obj.data.q + ...
                    randn(size(obj.data,1), 1) * obj.model.varepsilon;
            end
            mr = obj.means({'p', 'q'}, 'productid') ;
            if obj.model.endog
                obj.demand.var.instruments = sprintf('inst%d ', 1:6);
            end
            obj.demand.data = obj.data;
            obj.simDemand.data = obj.data;
        end
        
        function mr = simulateDemand(obj)
            obj.market.demand = obj.simDemand;
            if isempty(obj.model.firm)
                obj.market.var.firm = 'productid'; % One product firms
            else
                obj.market.var.firm = 'firm'; % Firm variable has been created
            end
            % Existence of var.costs name could be alternative to:
            obj.market.c = obj.data.c;
            
            obj.market.equilibrium();
            
            obj.data.p = obj.market.p;
            if obj.model.varepsilon > 0
                obj.simDemand.d = obj.simDemand.d + ...
                    randn(size(obj.data,1), 1) * obj.model.varepsilon;
                obj.simDemand.init();
            end
            obj.data.q = obj.simDemand.getDemand(obj.data.p);
           
            mr = obj.means({'p', 'q'}, 'productid') ;
            if obj.model.endog && isempty(obj.demand.var.instruments)
                if obj.model.randomProducts
                    obj.demand.var.instruments = 'nprod';
                elseif obj.model.gamma ~= 0
                    obj.demand.var.instruments = 'w';
                else
                    obj.demand.var.instruments = 'c';
                end
            end
            obj.demand.data = obj.data;
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
        function randDraws(varargin)
            if nargin > 0
                seed = varargin{1};
            else
                seed = 9;
            end
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed', seed));
        end
        
        function vals = testSame(x, y, varargin)
            if nargin > 2
                assert(all( abs(x - y) < 10e-6), 'Test %d failed' , varargin{1});
            else
                assert(all( abs(x - y) < 10e-6), 'Test failed');
            end
        end
        
        function vals = testEqual(x, y, sensitivity)
        % testEqual(x, y, sens) tests whether x and y are within 
        % sensitivity sens in percentage terms
            diff = abs(max((y - x) ./ y));
            if length(x) == 1
                vals = [y, x, diff];
            else
                vals = diff;
            end
            if diff > sensitivity
                disp(vals)
                error('True, calculated and percentage diff are greater than %f',...
                    sensitivity);
            end
        end
    end
end

