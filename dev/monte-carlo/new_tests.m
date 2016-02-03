% Put new test here while testing
% clear

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2);

m = SimMarket();
m.model.endog = false;
m.model.beta = [-1; 1; 2];

m.model.randproducts = false;
m.model.optimalIV = false;

%% Test: Optimal IV
m.randdraws();
m.demand = NestedLogitDemand;
% m.demand.var.nonlinear = 'x';

m.model.optimalIV = true;
m.model.endog = false;
m.model.randproducts = false;
% m.model.products = 10;
m.init();

results = m.simulateDemand()
m.simDemand.data = m.data;
m.simDemand.var.depvar = 'sh';
m.simDemand.init();
m.simDemand.X*m.simDemand.beta

return
results = m.estimate()
testtrue1(results{'p','Coef'} , m.model.beta(1) )

