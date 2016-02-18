%% Test 9: Testing ordering of multiple nonlinear parameters
display '**********************  Test 1  *************************'
for tst = 1:1
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true; 
m.demand.config.compiled = false;
% m.demand.config.tolerance = 1e-12;
m.demand.settings.drawmethod = 'quadrature';
m.demand.config.guessdelta = false;
m.demand.config.fpmaxit = 3000; 
m.demand.var.exog = 'x1 x2';
m.demand.alpha = 1;
m.demand.sigma = [.1; .2];
if ~logical(tst)
    m.demand.var.nonlinear = {'x1 p', 'normal'};
else
    m.demand.var.nonlinear = {'p x1', 'normal'};
end
m.demand.settings.nind = 1000;
%  m.demand.settings.marketdraws = true;
% m.demand.var.instruments = 'nprod nprod2 c';
% m.demand.settings.optimalIV = true;
% m.model.gamma = 1;
% m.model.endog = true;
% m.model.randomProducts = true;
% Increase in number of observations to get significance
m.model.markets = 200;
% m.model.products = 10;
m.model.x = [5, 0, 0];
m.model.x_vcv = [1, 1, 1];
m.model.beta = [ 1, -1, 20];

m.create();

display(m.model)

results = m.estimate()
sp = m.estDemand.results.sigma0
% SimMarket.testEqual(results{'rc_lP','Coef'} , results{'rc_lP', 'True_val'}, 1e-1 )

m.findCosts()
meanCosts = mean(m.data.c)
% SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )
end
return

 %% Test 9: Testing ordering of multiple nonlinear parameters
display '**********************  Test 2  *************************'
for tst = 0:1
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true; 
m.demand.config.compiled = false;
% m.demand.config.tolerance = 1e-12;
% m.demand.settings.drawmethod = 'quadrature';
m.demand.config.guessdelta = false;
m.demand.config.fpmaxit = 3000; 
m.demand.var.exog = 'x1 x2';
m.demand.alpha = 1;
if ~logical(tst)
    m.demand.sigma = [1; .2];
    m.demand.var.nonlinear = {'x1 x2', 'normal'};
else
    m.demand.sigma = [.2; 1];
    m.demand.var.nonlinear = {'x2 x1', 'normal'};
end
m.demand.settings.nind = 1000;
%  m.demand.settings.marketdraws = true;
% m.demand.var.instruments = 'nprod nprod2 c';
% m.demand.settings.optimalIV = true;
% m.model.gamma = 1;
% m.model.endog = true;
% m.model.randomProducts = true;
% Increase in number of observations to get significance
m.model.markets = 500;
% m.model.products = 10;
m.model.x = [5, 0, 0];
m.model.x_vcv = [1, 1, 1];
m.model.beta = [ 1, -1, 0];

m.create();

display(m.model)

results = m.estimate()
% SimMarket.testEqual(results{'rc_lP','Coef'} , results{'rc_lP', 'True_val'}, 1e-1 )

m.findCosts()
meanCosts = mean(m.data.c)
% SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )
end
 