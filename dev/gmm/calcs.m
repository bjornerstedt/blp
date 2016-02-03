% This calculation shows that the variance of costs does not depend much on
% beta. Variances are calculated at starting point beta and actual beta

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2);

m = SimMarket();
m.model.endog = false;
m.model.beta = [-1; 1; 2];

m.model.randproducts = false;
m.model.optimalIV = false;

m.demand = NestedLogitDemand;
% m.demand.var.nonlinear = 'x';

m.model.optimalIV = true;
m.model.endog = false;
m.model.randproducts = false;
% m.model.products = 10;
m.init();
m.simDemand.settings.paneltype = 'none';

results = m.simulateDemand()

% Now change the beta and calculate costs and variances
m.simDemand.beta = [-3; 1; 3];
m.simDemand.data = m.data;
m.simDemand.var.depvar = 'sh';
m.simDemand.init();
epsilon = m.simDemand.y - m.simDemand.X*m.simDemand.beta;
m.findCosts(m.simDemand)
avcosts = mean(m.market.c)
varcalc = [var(epsilon), var(m.market.c)]

% Calculate variances for estimate
results = m.estimate()
m.findCosts(m.estDemand)
mean(m.market.c)
epsilon = m.estDemand.y - m.estDemand.X*m.estDemand.beta;
varcalc = [var(epsilon), var(m.market.c)]



