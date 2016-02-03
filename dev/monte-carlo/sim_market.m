%% MONTE-CARLO NESTED AND MIXED LOGIT MODELS 
% Simulate data based on model and estimate
% Products exist with exog probability and count instruments based on the
% number of products are used to handle endogeneity

clear 

% Define model and est unless variable modelDef == true
if ~exist('modelDef', 'var') || ~modelDef

    % Note that both RC and instruments can be separated into model and
% estimation strategies by having  different vars in model and est.
model.RC = true;            % RC or NL
model.endog = true;        % Endog with count instruments or no endog
model.markets = 200;  
model.products = 5;
model.randproducts = true; % Let the number of products be random

model.beta = [-1; 1; 0];
model.rc_sigma = .1;
model.x = [5,0];
model.x_sigma = [1,1];
model.c = 4;
model.c_sigma = 1;

model.epsilon_sigma = .1;
model.sigma_xi = .1;
model.endog_sigma = 0.1;           % Degree of correlation between x and epsilon
model.prob_prod = .8;       % Probability of product existing in market

est.panel = true;           % LSDV panel or pooled estimate
est.method = 'hypercube';
est.nonlin = 'x';
est.optimalIV = true;
est.findCosts = false;          % Calculate costs or Simulate based on set costs
est.sameDraws = false;
end


%% Creation of demand data

rseed = @()RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 999));
rseed();
data = table();
n = model.markets * model.products;  % number of observations  
data.marketid = reshape(repmat(1:model.markets, model.products, 1), n, 1);
data.productid = repmat((1:model.products)', model.markets, 1);

% epsilon_jt = varepsilon_jt + xi_j
epsilon =  randn(n, 1) * model.epsilon_sigma + ...
    repmat(randn(model.products, 1) * model.sigma_xi, model.markets, 1);

if model.endog
    data.p = model.x(1) + randn(n, 1)*model.x_sigma(1) + model.endog_sigma*epsilon;
else
    data.p = model.x(1) + randn(n, 1)*model.x_sigma(1); 
end
data.x = model.x(2) + data.productid/model.products + randn(n, 1)*model.x_sigma(2);    
data.constant = ones(n, 1);
x0 = [table2array(data(:, 'x')), data.constant];
data.d = x0*model.beta(2:end) + epsilon;
data.c = data.constant*model.c + randn(n, 1)*model.c_sigma;

%% Random selection of products

if model.randproducts
    % probability prob_prod of a product existing in a period
    prodsel = logical(binornd(1, model.prob_prod, n, 1));
    data = data(prodsel, :);
end
% Create count instrument
nprod = accumarray(data.marketid, data.constant);
data.nprod = nprod(data.marketid,:);
data.nprod2 = data.nprod .^ 2;
%% Simulation of dataset
if model.RC
    demand = MixedLogitDemand(data);
    demand.settings.drawmethod = est.method;
    demand.var.nonlinear = est.nonlin;
    demand.beta = model.beta;
    demand.rc_sigma = model.rc_sigma;
else
    demand = NestedLogitDemand(data);
end
demand.var.market = 'marketid';
demand.var.price = 'p';
demand.var.exog = 'x';
demand.settings.paneltype = 'none'; % Makes no difference in simulation
demand.alpha = -model.beta(1);
demand.d = data.d;

if est.findCosts
    data.sh = zeros(size(data.d));
    for t = 1:model.markets
        selection = data.marketid == t;
        demand.initSimulation(t);
        data.sh(selection) = demand.shares(data.p(selection));
    end
    demand.share.s = data.sh;
    shtab = table(accumarray(data.productid, data.sh,[],@mean), ...
        'VariableNames', {'Shares'})
    display 'Average sum shares'
    disp(mean(accumarray(data.marketid, data.sh)))
end

market = Market(demand);
market.var.firm = 'productid'; % One product firms
market.c = data.c;
if est.findCosts
    market.findCosts();
    ctab = table(accumarray(data.productid, market.c,[],@mean), ...
        'VariableNames', {'Costs'})
    display 'Mean costs'
    mean(table2array(ctab))
else
    market.equilibrium();
    ptab = table(accumarray(data.productid, market.p,[],@mean), ...
        'VariableNames', {'Sim_price'})
    data.p = market.p;
    data.sh = market.s;
    
    %% Estimation
    if est.sameDraws
        rseed();
    end
    
    if model.RC
        demand = MixedLogitDemand(data);
        demand.settings.drawmethod = est.method;
        demand.var.nonlinear = est.nonlin;
    else
        demand = NestedLogitDemand(data);
    end
    demand.var.market = 'marketid';
    demand.var.price = 'p';
    demand.var.exog = 'x';
    demand.var.quantity = 'sh';
    demand.var.marketsize = 'constant';
    demand.settings.paneltype = 'none';
    if model.endog
        demand.var.instruments = 'nprod nprod2';
    end
    if est.panel
        demand.var.market = 'marketid';
        demand.var.panel = 'productid';
        demand.settings.paneltype = 'lsdv';
    else
        demand.settings.paneltype = 'none';
    end
    result = demand.estimate();
    if model.RC
        if model.endog && est.optimalIV
            demand.settings.optimalIV = true;
            result = demand.estimate();
        end
        truevals = table([model.beta; model.rc_sigma]);
        truevals.Properties.VariableNames = {'Theta'};
    else
        truevals = table(model.beta);
        truevals.Properties.VariableNames = {'Beta'};
    end
    
    disp([truevals, result]);
end
