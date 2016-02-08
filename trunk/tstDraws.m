clear 
% Simulated Testing

testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);
%% Test 9: Exogenous nonlinear price constant, testing findCosts
display '**********************  Test 9  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.var.nonlinear = 'p x';
m.demand.alpha = 1;
m.demand.sigma = [0.5; 1];
m.demand.settings.nind = 500;

% m.demand.settings.optimalIV = false;
% m.model.endog = false;

% Increase in number of observations to get significance
m.model.markets = 200;
m.model.products = 10;
m.create();
display(m.model)

results = m.estimate()
SimMarket.testEqual(results{'rc_p','Coef'} , results{'rc_p', 'True_val'}, 3e-2 )

m.findCosts()
meanCosts = mean(m.data.c)
SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )
