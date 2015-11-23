% Creation of demand data
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 999));
            
% m = SimMarket(NestedLogitDemand);
m = SimMarket(MixedLogitDemand);
m.model.endog = true;
m.model.beta = [-.2; 1; -5];
m.model.x = [5,5];
m.model.markets = 100;  
m.model.randproducts = true; 
% m.model.products = 5;


m.init();
if true
    m.calculateDemand()
    m.findCosts(m.simDemand)
else
    m.simulateDemand()
end
m.estimate()
m.findCosts(m.estDemand)
