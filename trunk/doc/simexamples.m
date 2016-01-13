clear
%% NL Demand Simulation
% * Only estimate

demand = NLDemand();
demand.alpha = 1;
m1 = SimMarket()
m1.demand = demand;

m1.model.markets = 100;

m1.create()
dt1 = m1.data;

%% NL Demand estimation

demand = NLDemand(dt1);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';
demand.estimate()


%% RC Demand simulation
% * random products
% * multi product firm
% * Merger

demand2 = RCDemand();
demand2.alpha = .3;
demand2.sigma = 1;
demand2.var.nonlinear = 'x';

m2 = SimMarket();
m2.demand = demand2;

m2.model.endog = true;
m2.model.randomProducts = true;
m2.model.firm = [1,1,2,2,3];
m2.model.markets = 200;
m2.create();

dt2 = m2.data;
%% RC Demand estimation

demand2 = RCDemand(dt2);

demand2.var.market = 'marketid';
demand2.var.panel = 'productid';
demand2.var.price = 'p';
demand2.var.quantity = 'q';
demand2.var.marketsize = 'constant';
demand2.var.exog = 'x';
demand2.var.instruments = 'nprod nprod2 c';

demand2.var.nonlinear = 'x';
demand2.settings.optimalIV = true;

result = demand2.estimate()

market = Market(demand2);
market.var.firm = 'firm';

market.findCosts(dt2.marketid == 1);

% Note that summary is of one market
market.summary()

market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();

compare(market, market2)
compare(market, market2, 'selection', dt2.marketid == 1)

%% NL Demand with cost estimation

demand = NLDemand();
demand.alpha = 0.5;
demand.sigma = 0.5;
% demand.var.nests = 'type';

m3 = SimMarket();
m3.demand = demand;
m3.model.types = 2;

m3.market = Market;
m3.market.settings.conduct = 0.5;
m3.model.markets = 200;
m3.model.firm = [1,1,1,2,2];

m3.model.endog = true;
m3.model.randomProducts = true;

m3.model.gamma = 1;
m3.create();
dt3 = m3.data;

demand = NLDemand(dt3);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';

demand.var.nests = 'type';
demand.var.instruments = 'c w nprod';
result = demand.estimate()

market = Market(demand);
market.var.firm = 'firm';
market.findCosts();
market.y = market.c;
market.var.exog = 'w';
market.var.panel = 'productid';
costEstimate = market.estimate()