% Relevant tests for modification of MixedLogitDemandNew
clear

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2, ...
    'Percentage difference is %f', abs((x - y)/y));
testtrue5 = @(x,y)assert(abs((x - y)/y)<5e-2, ...
    'Percentage difference is %f', abs((x - y)/y));
demandArray = {MixedLogitDemand, MixedLogitDemand2};
demandSelect = 2;

m = SimMarket();
model = m.model;
model.endog = false;
model.markets = 100;
model.randproducts = false;
model.optimalIV = false;

% 
% %% Test: MixedLogitDemandNew
% m = SimMarket(model);
% m.demand = demandArray{demandSelect};
% m.init();
% results = m.simulateDemand()
% testdiff(results{1, 'sh'} , 0.016533)
% testdiff(results{1, 'p'} , 5.0648)
% results = m.estimate()
% testdiff(results{'p','Coef'} , -1.001)
% testdiff(results{'p','Std_err'} , 0.0050055)
% 
% testtrue(results{'p','Coef'} , m.model.beta(1) )

%% Test: Nonlinear price
m = SimMarket(model);
m.demand = demandArray{demandSelect};
m.model.nonlinear = 'x';
m.model.rc_sigma = [ 1];

m.model.optimalIV = false;
m.model.endog = true;
m.model.randproducts = true;
m.model.markets = 100;
m.model.products = 10;
m.model.x_sigma = [.2, 1];
m.model.c = 0.9;
m.init();

results = m.simulateDemand()
m.data.nprod3 = m.data.nprod.*m.data.x;
m.data.nprod4 = m.data.x.^2;
m.data.nprod5 = m.data.nprod2.*m.data.x.^2;
m.estDemand.config.hessian = true;
m.simDemand.data = m.data;
m.simDemand.var.instruments = 'nprod nprod2 nprod3 nprod4 nprod5'
results = m.estimate()
testtrue5(results{'p','Coef'} , m.model.beta(1) )