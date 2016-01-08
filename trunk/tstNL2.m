%% Test 2: NestedLogitDemand - two level nests
display '**********************  Test 2  *************************'
m = SimMarket();
m.demand = NestedLogitDemand;
m.demand.var.nests = 'type1 type2';
m.demand.alpha = 2;
m.demand.sigma = [0.6, 0.3];
m.model.types = [2, 3];
m.model.products = 9;
m.create();

display(m.model)

results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'}, 1e-2)

testVals = [results{'p','Coef'}, results{'lsjh','Coef'}, results{'p','Std_err'}];
correctVals = [-1.995650385534655,0.502068400039377,0.006222691020318];

testSameResults(testVals, correctVals, 2);