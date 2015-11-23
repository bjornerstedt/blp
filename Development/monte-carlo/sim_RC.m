%% RANDOM COEFFICIENT LOGIT MODEL WITH SIMULATED DATA 

clear 
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 999));

instruments = true
optimalIV = true
markets = 100;  
products = 20;
panel = true;
sigma_epsilon = .1;
sigma_xi = .1;
beta = [-.2; 1; -5];
rc_sigma = 1;
method = 'hypercube'
nonlin = 'constant'
sameDraws = false

%% Creation of dataset

data = table();
n = markets*products;  % number of observations  
data.marketid = reshape(repmat(1:markets, products, 1), n, 1);
data.productid = repmat((1:products)', markets, 1);

epsilon = randn(n, 1) * sigma_epsilon;  
if panel
    % Random effects version unless instruments=true
    xi = repmat(randn(products, 1)*sigma_xi, markets, 1);
    epsilon = epsilon + xi;
end

% Simulate Instruments a la Nevo
if instruments
    A = 2 + 0.2 * randn(n, 6);
    M = ones(6, 6)*0.8;
    M(1:7:36) = 1; % Make diagonal elements = 1
    inst = A*chol(M);
    % Keep expected effect of instruments on p to be zero
    suminst = sum(inst, 2) - mean(sum(inst, 2)); 
    data.p = 5 + randn(n, 1)*2 + suminst + 0.5*epsilon; 
    data = [data, array2table(inst)];    
else
    data.p = 5 + randn(n, 1); 
end
data.x = 5 + randn(n, 1);    
data.constant = ones(n, 1);
x0 = [table2array(data(:, 'x')), data.constant];
d = x0*beta(2:end) + epsilon;

%% Simulation of dataset

demand = MixedLogitDemand(data);
demand.settings.drawmethod = method;
demand.var.nonlinear = nonlin;
demand.beta = beta;
demand.rc_sigma = rc_sigma; 
demand.var.market = 'marketid';
demand.var.price = 'p';
demand.var.exog = 'x';
demand.settings.paneltype = 'none'; % Makes no difference in simulation
demand.alpha = -beta(1);
demand.d = d;
demand.init();

data.sh = zeros(size(d));
for t = 1:markets
    selection = data.marketid == t;
    demand.initSimulation(t);
    data.sh(selection) = demand.shares(data.p(selection));
end
sh = data.sh(1:products)
sum(sh)

%% Estimation
if sameDraws
    RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 10));
end
demand = MixedLogitDemand(data);
demand.settings.drawmethod = method;
demand.var.nonlinear = nonlin;
demand.var.market = 'marketid';
demand.var.price = 'p';
demand.var.exog = 'x';
demand.var.quantity = 'sh';
demand.var.marketsize = 'constant';
demand.settings.paneltype = 'none';
if instruments
    demand.var.instruments = 'inst1 inst2 inst3 inst4 inst5 inst6';
end
if panel
    demand.var.market = 'marketid';
    demand.var.panel = 'productid';
    demand.settings.paneltype = 'lsdv';
else
    demand.settings.paneltype = 'none';
end
demand.init(); 

result = demand.estimate();
if instruments && optimalIV
    demand.settings.optimalIV = true;
    demand.estimate();
end
truevals = table([beta; rc_sigma]);
truevals.Properties.VariableNames = {'Theta'};
disp([truevals, result]);


