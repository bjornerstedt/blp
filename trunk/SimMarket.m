classdef SimMarket < matlab.mixin.Copyable 
% SimMarket creates simulated data and estimates on it
% MONTE-CARLO NESTED AND MIXED LOGIT MODELS 
% Simulate data based on model and estimate
% Products exist with exog probability and count instruments based on the
% number of products are used to handle endogeneity
%   $Id: SimMarket.m 142 2015-10-09 13:07:17Z d3687-mb $
    
    properties
        model
        data
        demand
        simDemand % Class with settings for data generation
        estDemand % Settings for estimation
        market
    end
    properties 
        means = @(T,x,ind)array2table(splitapply(@mean, T{:,x}, T.(ind)), ...
            'VariableNames',x);
        approx = @(x,y)(abs(x - y) < 10^-4)
        randdraws = @()RandStream.setGlobalStream(RandStream('mt19937ar','Seed', 999));

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
                obj.model.endog = false;      % Endog with count instruments or no endog
                obj.model.markets = 200;
                obj.model.products = 5;
                obj.model.types = 1;
                obj.model.randproducts = true; % Let the number of products be random
                
                % p and x mean and sd:
                obj.model.x = [5,0];
                obj.model.x_sigma = [1,1];
                obj.model.c = 4;
                obj.model.c_sigma = 1;
                obj.model.gamma = 0;
                
                % Individual and product level shocks:
                obj.model.epsilon_sigma = .1;
                obj.model.sigma_xi = .1;
                
                obj.model.endog_sigma = 0.1; % Degree of corr between x and epsilon
                obj.model.prob_prod = .8;    % Prob of product existing in market
                
                obj.model.alpha = -1;
                obj.model.sigma = [];
                obj.model.rc_sigma = 1;
                obj.model.beta = [ 1; 0];
                
                obj.model.optimalIV = true;
                obj.model.nonlinear = 'constant';
                obj.model.drawmethod = 'hypercube';
                obj.model.nests = '';
            end
        end
        
        % Separate from constructor to allow init from model 
        function init(obj)
            % Simulation of dataset
            if isa(obj.demand, 'MixedLogitDemand')
                if ~isempty(obj.model.drawmethod)
                    obj.demand.settings.drawmethod = obj.model.drawmethod;
                end
                if ~isempty(obj.model.nonlinear)
                    obj.demand.var.nonlinear = obj.model.nonlinear;
                end
                obj.demand.rc_sigma = obj.model.rc_sigma;
            end
            if isempty(obj.demand.var.nests)
                obj.demand.var.nests = obj.model.nests;
            end
            
            obj.demand.var.market = 'marketid';
            obj.demand.var.panel = 'productid';
            obj.demand.var.price = 'p';
            obj.demand.var.exog = 'x';
            obj.demand.var.quantity = 'q';
            obj.demand.var.marketsize = 'constant';
            % Not used in simulation
            obj.demand.settings.paneltype = 'lsdv';
            
            obj.simDemand = copy(obj.demand);
            obj.estDemand = copy(obj.demand);

            obj.createData();
            obj.simDemand.data = obj.data;
            obj.simDemand.d = obj.data.d; % <<<<<<<<<<< init is required
            
            obj.simDemand.beta = [obj.model.alpha; obj.model.beta];            
            obj.simDemand.alpha = -obj.model.alpha;
            if isa(obj.demand, 'NestedLogitDemand')
                obj.simDemand.sigma = obj.model.sigma;
            end
        end
        
        function createData(obj)
            obj.data = table();
            n = obj.model.markets * obj.model.products;    
            obj.data.marketid = reshape(repmat(1:obj.model.markets, ...
                obj.model.products, 1), n, 1);
            obj.data.productid = repmat((1:obj.model.products)', ...
                obj.model.markets, 1);
            if obj.model.types > 1
                reps = ceil(obj.model.products / obj.model.types);
                typelist = repmat((1:obj.model.types)', reps, 1);
                obj.data.type = repmat(typelist(1:obj.model.products), ...
                    obj.model.markets, 1);
            end

            % epsilon_jt = varepsilon_jt + xi_j
            epsilon =  randn(n, 1) * obj.model.epsilon_sigma + ...
                repmat(randn(obj.model.products, 1) * obj.model.sigma_xi, ...
                obj.model.markets, 1);
 
            % Create random price:
            if obj.model.endog
                % Instruments a la Nevo:
                ninstr = 6; 
                A = 2 + 0.2 * randn(n, ninstr);
                M = eye(ninstr)*0.2 + ones(ninstr, ninstr)*0.8;
                inst = A * chol(M);
                % Keep expected effect of instruments on p to be zero
                suminst = sum(inst, 2) - mean(sum(inst, 2));
                obj.data.p = obj.model.x(1) + randn(n, 1) * obj.model.x_sigma(1) ...
                    + suminst + obj.model.endog_sigma * epsilon;
            else
                obj.data.p = obj.model.x(1) + randn(n, 1) * obj.model.x_sigma(1);
                inst = [];
            end
           
            % Create random x var:
            obj.data.x = obj.model.x(2) + obj.data.productid / obj.model.products ...
                + randn(n, 1) * obj.model.x_sigma(2);    
            obj.data.constant = ones(n, 1);
            
            % This should be done in the demand classes:            <<<<<<<<<<<<<<<<
            x0 = [table2array(obj.data(:, 'x')), obj.data.constant];
            obj.data.d = x0 * obj.model.beta + epsilon;
            
            if obj.model.gamma ~= 0 % Otherwise existing tests fail
                obj.data.w = randn(n, 1);
            else
                obj.data.w = zeros(n,1);
            end
            obj.data.c = obj.data.constant * obj.model.c ...
                + obj.model.gamma * obj.data.w ...
                + randn(n, 1) * obj.model.c_sigma;

            obj.data = [obj.data, array2table(inst)];
            
            % Random selection of products
            if obj.model.randproducts
                % probability prob_prod of a product existing in a period
                prodsel = logical(binornd(1, obj.model.prob_prod, n, 1));
                obj.data = obj.data(prodsel, :);
            end
            % Create count instrument
            nprod = accumarray(obj.data.marketid, obj.data.constant);
            obj.data.nprod = nprod(obj.data.marketid,:);
            obj.data.nprod2 = obj.data.nprod .^ 2;  
            % Add instruments last
        end
        
        % Function can reset data and i
        function mr = calculateDemand(obj)
            obj.simDemand.init();
            obj.data.q = obj.simDemand.getDemand(obj.data.p);
            
            mr = obj.means(obj.data, {'p', 'q'}, 'productid') ;
            
            display 'Average sum shares'
            disp(mean(accumarray(obj.data.marketid, obj.data.q)))
            if obj.model.endog
                obj.estDemand.var.instruments = sprintf('inst%d ', 1:6);
            end
            obj.estDemand.data = obj.data;
            obj.simDemand.data = obj.data;
        end
        
        function mr = findCosts(obj, demand)            
            obj.market = Market(demand);
            obj.market.var.firm = 'productid'; % One product firms
            obj.market.findCosts();
            mr = obj.means(obj.data, {'sh','c',  'p'}, 'productid') ;
            cr = rowfun(@(x,y)(y-x)/y, mr(:,2:3),'OutputVariableNames','Markup');
            mr = [mr,cr];
        end
        
        function mr = simulateDemand(obj, varargin)
            obj.simDemand.init();
            if nargin > 1
                obj.market = Market(varargin{1});
            else
                obj.market = Market(obj.simDemand);                
            end
            obj.market.var.firm = 'productid'; % One product firms
            % Existence of var.costs name could be alternative to:
            obj.market.c = obj.data.c;

            obj.market.equilibrium();
            
            obj.data.p = obj.market.p;
            obj.data.sh = obj.market.s;
            obj.data.q = obj.simDemand.getDemand(obj.data.p);    
            
            mr = obj.means(obj.data, {'p', 'sh'}, 'productid') ;
            if obj.model.endog
                if obj.model.randproducts
                    obj.estDemand.var.instruments = 'nprod nprod2';
                else
                    obj.estDemand.var.instruments = 'c';
                end
            end
            obj.estDemand.data = obj.data;
       end
                   
       function results = estimate(obj)           
            result = obj.estDemand.estimate();
            if isa(obj.demand, 'NestedLogitDemand')
                beta = [obj.model.alpha; obj.model.sigma; obj.model.beta];
            else
                beta = [obj.model.alpha; obj.model.beta];
            end
            if strcmpi(obj.estDemand.settings.paneltype, 'fe')
                beta = beta(1:end-1);
            end
            if isa(obj.demand, 'MixedLogitDemand')
                if obj.model.endog && obj.model.optimalIV
                    obj.estDemand.settings.optimalIV = true;
                    result = obj.estDemand.estimate();
                end
                truevals = table([beta; obj.model.rc_sigma]);
                truevals.Properties.VariableNames = {'Theta'};
            else
                truevals = table(beta);
                truevals.Properties.VariableNames = {'Beta'};
            end
            results = [truevals, result];
        end
        
        function results = run(obj)
            obj.createData();
            obj.initDemand();
            obj.simulateDemand();
            if obj.est.findCosts
                obj.findCosts();
            else
                results = obj.estimate();
            end
        end
    end
    
end

