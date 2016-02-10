clear 
% Simulated Testing

testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);

%% Test 7: Optimal IV
display '**********************  Test 7  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.var.nonlinear = 'x';
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.settings.optimalIV = true;

m.model.endog = true;
m.model.randomProducts = true;
m.model.products = 10;
m.create();
display(m.model)

results = m.estimate()

testVals = [results{'p','Coef'}, results{'rc_x','Coef'}];
correctVals = [-0.969984767147493,1.044593651632106];
testSameResults(testVals, correctVals, 7);

SimMarket.testEqual(results{'rc_x','Coef'} , results{'rc_x', 'True_val'}, 1e-1 )
