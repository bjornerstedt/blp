clear
%% Test 1: NLDemand with findCosts 
% Test that findCosts and estimate are inverses
display '**********************  Test 1  *************************'

m1 = SimMarket()
demand = NLDemand();
demand.alpha = .3;
m1.demand = demand;

m1.model.endog = true;
m1.model.markets = 100;
m1.model.products = 3;
m1.model.beta = [ 1, 0];

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
% demand.var.instruments = 'nprod nprod2';
demand.var.instruments = 'c';
% demand.var.instruments = 'inst1 inst2 inst3 inst4 inst5 inst6';
demand.estimate()

market = Market(demand);
market.var.firm = 'productid';

display 'Test that mean calculated costs are close to actual'
market.findCosts( );
cc = [dt1.c, market.c];
% meanCosts = SimMarket.testEqual(mean(dt1.c), mean(market.c), 1)

display 'Test that mean equilibrium prices are close to starting vals'
market.equilibrium( )
pp = [dt1.p, market.p];
% meanPrices = SimMarket.testEqual(mean(dt1.p), mean(market.p), 1)

market.findCosts( );
cc = [cc, market.c];
market.equilibrium( );
pp = [pp, market.p];
market.summary()
market.summary('selection', dt1.marketid == 1)
market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();

compare(market, market2)
compare(market, market2, 'selection', dt1.marketid == 1)

%% RCDemand simulation

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

result = demand2.estimate()

market = Market(demand2);
market.var.firm = 'firm';

market.findCosts();
averageCosts = mean(market2.c)
market.summary()
market.summary('selection', dt2.marketid == 1)

market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();

compare(market, market2)
compare(market, market2, 'selection', dt2.marketid == 1)
