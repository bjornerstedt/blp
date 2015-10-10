classdef SimMarket  < matlab.mixin.Copyable 
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
                obj.model.RC = true;          % RC or NL
                obj.model.endog = false;      % Endog with count instruments or no endog
                obj.model.ces = false;
                obj.model.markets = 200;
                obj.model.products = 5;
                obj.model.randproducts = true; % Let the number of products be random
                
                obj.model.beta = [-1; 1; 0];
                obj.model.rc_sigma = 1;
                
                obj.model.x = [5,0];
                obj.model.x_sigma = [1,1];
                obj.model.c = 4;
                obj.model.c_sigma = 1;
                
                obj.model.epsilon_sigma = .1;
                obj.model.sigma_xi = .1;
                obj.model.endog_sigma = 0.1; % Degree of corr between x and epsilon
                obj.model.prob_prod = .8;    % Prob of product existing in market
                
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
            obj.simDemand.d = obj.data.d;
            if obj.model.ces
                obj.simDemand.settings.ces = true;
                obj.estDemand.settings.ces = true;
            end
            
            obj.simDemand.beta = obj.model.beta;            
            obj.simDemand.alpha = -obj.model.beta(1);
            
        end
        
        function createData(obj)
            obj.data = table();
            n = obj.model.markets * obj.model.products;    
            obj.data.marketid = reshape(repmat(1:obj.model.markets, ...
                obj.model.products, 1), n, 1);
            obj.data.productid = repmat((1:obj.model.products)', ...
                obj.model.markets, 1);

            % epsilon_jt = varepsilon_jt + xi_j
            epsilon =  randn(n, 1) * obj.model.epsilon_sigma + ...
                repmat(randn(obj.model.products, 1) * obj.model.sigma_xi, ...
                obj.model.markets, 1);

            if obj.model.endog
                A = 2 + 0.2 * randn(n, 6);
                M = ones(6, 6)*0.8;
                M(1:7:36) = 1; % Make diagonal elements = 1
                inst = A*chol(M);
                % Keep expected effect of instruments on p to be zero
                suminst = sum(inst, 2) - mean(sum(inst, 2));
                obj.data.p = obj.model.x(1) + randn(n, 1)*obj.model.x_sigma(1) ...
                    + suminst + obj.model.endog_sigma * epsilon;
            else
                obj.data.p = obj.model.x(1) + randn(n, 1)*obj.model.x_sigma(1);
                inst = [];
            end
            
            obj.data.x = obj.model.x(2) + obj.data.productid/obj.model.products ...
                + randn(n, 1)*obj.model.x_sigma(2);    
            obj.data.constant = ones(n, 1);
            x0 = [table2array(obj.data(:, 'x')), obj.data.constant];
            obj.data.d = x0*obj.model.beta(2:end) + epsilon;
            obj.data.c = obj.data.constant * obj.model.c ...
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
            obj.data.sh = zeros(size(obj.data.d));
            obj.simDemand.init();
            % This functionality should be in NestedLogitDemand:
            for t = 1:obj.model.markets
                selection = obj.data.marketid == t;
                obj.simDemand.initSimulation(t);
                obj.data.sh(selection) = obj.simDemand.shares(obj.data.p(selection));
            end            
            if obj.model.ces
                obj.data.q = obj.data.sh ./ obj.data.p;
            else
                obj.data.q = obj.data.sh;
            end
            mr = obj.means(obj.data, {'p', 'sh'}, 'productid') ;
            
            display 'Average sum shares'
            disp(mean(accumarray(obj.data.marketid, obj.data.sh)))
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
            if obj.model.ces
                obj.data.q = obj.data.sh./obj.data.p;
            else
                obj.data.q = obj.data.sh;
            end
            mr = obj.means(obj.data, {'p', 'sh'}, 'productid') ;
            if obj.model.endog
                obj.estDemand.var.instruments = 'nprod nprod2';
            end
            obj.estDemand.data = obj.data;
       end
        
        function results = estimate(obj)           
            result = obj.estDemand.estimate();
            if strcmpi(obj.estDemand.settings.paneltype, 'fe')
                beta = obj.model.beta(1:end-1);
            else
                beta = obj.model.beta;
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

