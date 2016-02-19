%% Test 1: Testing ordering of multiple nonlinear with price
display '**********************  Test 1  *************************'
endog = true;
for tst = 0:1
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true; 
m.demand.config.compiled = false;
m.demand.settings.drawmethod = 'quadrature';
m.demand.config.guessdelta = false;
m.demand.config.fpmaxit = 5000; 
m.demand.alpha = .8;
if ~logical(tst)
    m.demand.sigma = [1; .2];
    m.demand.var.nonlinear = 'x p';
else
    m.demand.sigma = [.2; 1];
    m.demand.var.nonlinear = 'p x';
end
m.demand.settings.nind = 1000;
if endog
    m.demand.var.instruments = 'nprod nprod2 c';
    m.demand.settings.optimalIV = true;
    m.model.endog = true;
    m.model.randomProducts = true;
end
m.model.markets = 500;
% m.model.products = 10;
% m.model.x = [5, 0, 0];
% m.model.x_vcv = [1, 1, 1];
% m.model.beta = [ 1, -1, 20];

m.create();

display(m.model)

results = m.demand.estimate()
sp = m.demand.results.sigma0
% SimMarket.testEqual(results{'rc_lP','Coef'} , results{'rc_lP', 'True_val'}, 1e-1 )

m.findCosts()
meanCosts = mean(m.data.c)
% SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )
end
return

 %% Test 2: Testing ordering of multiple nonlinear parameters
display '**********************  Test 2  *************************'
for tst = 0:1
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true; 
m.demand.config.compiled = false;
% m.demand.config.tolerance = 1e-12;
m.demand.settings.drawmethod = 'quadrature';
% m.demand.config.guessdelta = false;
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

results = m.demand.estimate()
% SimMarket.testEqual(results{'rc_lP','Coef'} , results{'rc_lP', 'True_val'}, 1e-1 )
m.demand.results.other.fval
m.findCosts()
meanCosts = mean(m.data.c)
% SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )
end
 
%% Test 3: CES RCDemand
display '**********************  Test 3  *************************'
m = SimMarket();
m.demand = RCDemand;
% m.demand.settings.ces = true;
%         m.demand.settings.robust = 'false';
m.demand.alpha = 4;
m.demand.sigma = .5;
m.demand.settings.sigma0 = [1];
m.model.beta = [ 1, 4];
% m.model.markets = 200;
m.demand.var.nonlinear = 'p';
% m.model.endog = true;
% m.model.randomProducts = false;
% m.model.pricesFromCosts = false;

m.create();
display(m.model)

results = m.demand.estimate()
sigma0 = m.demand.results.sigma0 

SimMarket.testEqual(results{'lP', 'Coef'} , results{'lP', 'True_val'}, 1e-2)

%% Test 4: Test setting sigma0
display '**********************  Test 4  *************************'
sigma0 = [];
for i = 1:2
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.drawmethod = 'quadrature';
m.demand.alpha = 1;
m.demand.sigma = [1];

m.demand.var.nonlinear = 'x';
if ~isempty(sigma0)
    m.demand.settings.sigma0 = sigma0;
end
m.demand.settings.paneltype = 'lsdv';
m.demand.var.instruments = 'nprod nprod2 c';
m.model.gamma = 1;
m.model.endog = true;
m.model.randomProducts = true;
m.model.markets = 200;
m.model.products = 5;
m.create();
display(m.model)
results = m.demand.estimate()
if i == 2
SimMarket.testEqual(results{'rc_x','Coef'} , firstsigma, 1e-7 )
end
sigma0 = m.demand.results.sigma0 
firstsigma = results{'rc_x','Coef'};
end
