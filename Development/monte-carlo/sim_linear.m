%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LINEAR MODEL WITH SIMULATED DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 999));

instruments = false;
sigma_epsilon = 1;
markets = 10;  
products = 5;
panel = true;
sigma_xi = 1;

%% SIMULATION OF DATASET

data = table();
n = markets*products;  % number of observations  
data.marketid = reshape(repmat(1:markets, products, 1), n, 1);
data.productid = repmat((1:products)', markets, 1);
beta = [ -1; 1; 10];

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
    % xi or epsilon?:
    data.p = 5 + randn(n, 1)*2 + sum(inst, 2)+ 0.5*epsilon; 
    data = [data, array2table(inst)];    
else
    data.p = 5 + randn(n, 1)*2; 
end
data.x = 5 + randn(n, 1)*2;    

constant=ones(n, 1);
X= [table2array( data(:, 3:end)), constant];
data.y = X * beta + epsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

demand = Estimate(data);
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
demand.var.depvar = 'y';
demand.var.exog = 'p x';
demand.init(); 

demand.estimate()



