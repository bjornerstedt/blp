clear
m1 = SimMarket()
demand = NLDemand();
demand.alpha = .1;
m1.demand = demand;
m1.model.endog = true;
m1.model.randproducts = false;
m1.create()
dt1 = m1.data;
return
demand = NLDemand(dt1);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';
demand.var.instruments = 'nprod nprod2';

demand.estimate()

market = Market(demand);
market.var.firm = 'productid';

market.findCosts();
cc = [dt1.c, market.c];
mean(cc)
market.equilibrium()
averageCosts = mean(market.c)
market.summary()
market.summary('selection', dt1.marketid == 1)
market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();

compare(market, market2)
compare(market, market2, 'selection', dt1.marketid == 1)
