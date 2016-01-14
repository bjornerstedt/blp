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
m.model.firm = [1,1,2,3,1,2,2,3]
m.create()

market = 'marketid';
typeNames = {'type1', 'type2'};
firm = {'firm'};
% firm = [];

% demand.data has to be set.
[instruments, names] = m.demand.countInstruments(market, [firm, typeNames], {'constant', 'type1'}, m.data)

