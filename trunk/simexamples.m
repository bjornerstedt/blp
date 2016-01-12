clear
%% Test 1: NLDemand with findCosts 
display '**********************  Test 1  *************************'

m1 = SimMarket()
demand = NLDemand();
demand.alpha =.1;
m1.demand = demand;
m1.model.simulatePrices = false;
m1.model.endog = false;
m1.model.randproducts = false;
% m1.model.gamma = 2;
% m1.model.x_vcv = [1, .01];
m1.model.markets = 100;
m1.model.products = 5;
m1.model.beta = [ 1, 2];

m1.create()
dt1 = m1.data;

demand = NLDemand(dt1);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';
% demand.var.instruments = 'nprod nprod2';
% demand.var.instruments = 'inst1 inst2 inst3 inst4 inst5 inst6';

demand.estimate()

market = Market(demand);
market.var.firm = 'productid';

display 'Test that mean calculated costs are close to actual'
market.findCosts();
% meanCosts = SimMarket.testEqual(mean(dt1.c), mean(market.c), 3e-2)

display 'Test that mean equilibrium prices are close to starting vals'
market.equilibrium()
meanPrices = SimMarket.testEqual(mean(dt1.p), mean(market.p), 1e-2)

% market.summary()
% market.summary('selection', dt1.marketid == 1)
% market2 = copy(market);
% market2.firm(market2.firm == 2 ) = 1;
% market2.equilibrium();
% 
% compare(market, market2)
% compare(market, market2, 'selection', dt1.marketid == 1)
