clear

m = SimMarket();
model = m.model;
model.markets = 25;

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
m.model.randproducts = false;
m.model.endog = true;
m.demand = NestedLogitDemand;
m.demand.settings.paneltype = 'none';
m.init();
m.simulateDemand();
m.estDemand.settings.paneltype = 'none';

results = m.estimate()

gmm1(m.estDemand)