clear

m = SimMarket();
model = m.model;
model.markets = 25;
model.epsilon_sigma = .1;
model.gamma = 1;

%% Test: Count instruments
m = SimMarket(model);
m.model.randproducts = true;
m.model.endog = true;
m.demand = NestedLogitDemand;
m.init();
m.simulateDemand();
results = m.estimate()

%% Test: NestedLogitDemand
m = SimMarket(model);
m.model.randproducts = true;
m.model.endog = true;
m.demand = NestedLogitDemand;
%m.demand.settings.paneltype = 'none';
m.init();
m.simulateDemand();

display('2SLS estimate')
m.estDemand.var.instruments = 'w';
results = m.estimate();

display('GMM estimate')
%m.estDemand.settings.paneltype = 'none';
m.estDemand.settings.estimateMethod = 'gmm';
results = m.estimate();

gmm_funcs(m.estDemand);
 