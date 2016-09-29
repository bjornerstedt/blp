% Count instruments: count product by date
% num, 
% numg, form
% numhg, substance form
% numf, firm
% numfg, firm form
% numfgh, firm substance form
% 
% missing:
% substance
% firm substance

m = SimMarket;
m.demand = NLDemand;
m.demand.alpha = 1;
m.model.markets = 5;
m.model.products = 8;
m.model.types = [2,3];
m.model.randomProducts = true;
m.model.firm = [1,1,2,3,1,2,2,3];
m.create()

market = 'marketid';
typeNames = {'type1', 'type2'};
firm = {'firm'};
% firm = [];

% instruments = Estimate.countInstruments(m.data, market, [firm, typeNames], {'constant', 'type1'});
instruments = Estimate.countInstruments(m.data, market, [firm, typeNames]);

demand = NLDemand();
demand.alpha = 0.5;
demand.sigma = 0.5;
demand.var.nests = 'type';

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
dt3.inst = Estimate.countInstruments(dt3, 'marketid', {'firm', 'type'});

demand = NLDemand(dt3);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketSize = 'constant';
demand.var.exog = 'x';

demand.var.nests = 'type';
demand.var.instruments = 'nprod';
result = demand.estimate()
