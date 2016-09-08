clear

m = SimMarket();
model = m.model;
model.markets = 25;
model.epsilon = 1;
model.gamma = 1;
% model.types = 2;
model.firm = [1,1,3,2,2]';

%% Test: Count instruments
m = SimMarket(model);
%m.model.randomProducts = true;
m.model.endog = true;
m.demand = NLDemand;
m.demand.alpha = 1;
m.create();
display Estimate
results = m.demand.estimate()

%% Test: NLDemand
% m = SimMarket(model);
% m.model.randomProducts = true;
% m.model.endog = true;
% m.demand = NLDemand;
% m.demand.alpha = 1;
% %m.demand.var.nests = 'type';
% 
% m.demand.settings.paneltype = 'fe';
% m.create()

display('2SLS estimate')
m.demand.var.instruments = 'w';
results = m.demand.estimate()


m.findCosts(m.demand)
mean(m.data.c)

display('GMM estimate')
%m.estDemand.settings.paneltype = 'none';
m.demand.settings.estimateMethod = 'gmm';
results = m.demand.estimate()

% Estimate alpha nonlinearly
% gmm_funcs(m.demand);

display('Simultaneous Estimate')
alpha = [-.3]';

market = Market(m.demand);
market.var.firm = 'productid';
market.settings.paneltype = 'none';
market.var.exog = 'w';
market.init();
    
theta = market.estimateGMM( alpha)

market.findCosts()
display 'Average market c:'
mean(market.c)
aa=market.summary()
sum(aa{:,5})