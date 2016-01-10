clear

m = SimMarket();
model = m.model;
model.markets = 25;
model.epsilon_sigma = .1;
model.gamma = 1;
model.types = 2;
%model.firm = [1,1,3,2,2]';

%% Test: Count instruments
m = SimMarket(model);
m.model.randproducts = true;
m.model.endog = true;
m.demand = NLDemand;
m.demand.alpha = 1;
m.create();
results = m.estimate()

%% Test: NLDemand
m = SimMarket(model);
m.model.randproducts = true;
m.model.endog = true;
m.demand = NLDemand;
m.demand.alpha = 1;
%m.demand.var.nests = 'type';

m.demand.settings.paneltype = 'fe';
m.create()

display('2SLS estimate')
m.estDemand.var.instruments = 'w';
results = m.estimate();

m.findCosts(m.estDemand)
mean(m.data.c)

display('GMM estimate')
%m.estDemand.settings.paneltype = 'none';
m.estDemand.settings.estimateMethod = 'gmm';
results = m.estimate();

% gmm_funcs(m.estDemand);

display('Simultaneous Estimate')
alpha = [-.3]';

market = Market(m.estDemand);
market.var.firm = 'productid';
market.settings.paneltype = 'none';
market.var.exog = 'w';
market.init();
    
beta = market.estimateGMM( alpha)
mean(market.c)

market.findCosts()
mean(market.c)
market.summary()