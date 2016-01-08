clear

% Put these last to avoid testi
testSameResults = @(x,y,i)true;
testtrue = @(x,y,z)true;

testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);
testtrue = @(x,y,z)assert(abs((x - y)/y)<z, 'Percentage diff is %f', abs((x - y)/y));

%% Test 1: NLDemand
display '**********************  Test 1  *************************'
m = SimMarket();
m.demand = NLDemand;
m.demand.alpha = 1;
sresults = m.create();
display(m.model)
results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'}, 1e-2)

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.022769720687512,4.895422144428684,-1.000390188323463,0.004938388080890];

testSameResults(testVals, correctVals, 1);

%% Test 2: NLDemand - one nest
display '**********************  Test 2  *************************'
m = SimMarket();
m.demand = NLDemand;
m.demand.var.nests = 'type';
m.demand.alpha = 2;
m.demand.sigma = 0.5;

m.model.types = 2;
m.create();
display(m.model)

results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'}, 1e-2)

testVals = [results{'p','Coef'}, results{'lsjg','Coef'}, results{'p','Std_err'}];
correctVals = [-1.995650385534655,0.502068400039377,0.006222691020318];

testSameResults(testVals, correctVals, 2);

%% Test 2b: NLDemand - two level nests
display '**********************  Test 2b *************************'
m = SimMarket();
m.demand = NLDemand;
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
correctVals = [-1.995280968011216,0.602548086049320,0.005338793317819];

testSameResults(testVals, correctVals, 2);

%% Test 3: RCDemand - rc_constant
% Estimate without instruments, with exogeneous p.
display '**********************  Test 3  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.var.nonlinear = 'constant';
sresults = m.create();
display(m.model)

results = m.estimate()
% Check that estimate different than initial guess:
assert(abs(m.estDemand.results.sigma0 - m.estDemand.sigma) > 0.5)

testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-2 )

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.027528441329171,4.922947217451795,-1.001769146394760,0.005015988697898];
testSameResults(testVals, correctVals, 3);

%% Test 4: RCDemand - rc_x
display '**********************  Test 4  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.var.nonlinear = 'x';
sresults = m.create();
display(m.model)

results = m.estimate()
testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.028966904126857,4.961760272683717,-1.001534413344814,0.005151635757488];
testSameResults(testVals, correctVals, 4);

%testtrue(results{'rc_x','Coef'} , m.demand.sigma, 1e-1 )

%% Test 5: CES RCDemand
display '**********************  Test 5  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true;
%         m.estDemand.settings.robust = 'false';
m.model.endog = false;
m.demand.alpha = 4;
m.demand.sigma = 1;
m.model.beta = [ 1; 4];
m.model.markets = 200;
m.model.randproducts = false;
m.demand.var.nonlinear = 'constant';
m.model.simulatePrices = false;

m.create();
display(m.model)

results = m.estimate()
testVals = [results{'lP', 'Coef'}, results{'rc_constant', 'Coef'}];
correctVals = [-4.005833090434905,0.980403597281313];
testSameResults(testVals, correctVals, 5);

testtrue(results{'lP', 'Coef'} , results{'lP', 'Theta'}, 2e-2)

%% Test 6: Optimal IV
display '**********************  Test 6  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.var.nonlinear = 'x';
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.settings.optimalIV = true;

m.model.endog = true;
m.model.randproducts = true;
m.model.products = 10;
m.create();
display(m.model)

results = m.estimate()

testVals = [results{'p','Coef'}, results{'rc_x','Coef'}];
correctVals = [-0.969984767147493,1.044593651632106];
testSameResults(testVals, correctVals, 6);

testtrue(results{'rc_x','Coef'} , results{'rc_x', 'Theta'}, 1e-1 )

%% Test 7: LSDV/FE
% Four tests, comparing LSDV/FE for NL/FE
display '**********************  Test 7  *************************'
results = cell(2,1);
for d = 1:2
    paneltype = {'lsdv', 'fe'};
    for i = 1:2
        m = SimMarket();
       if d == 1
            m.demand = NLDemand;
        else
            m.demand = RCDemand;
            m.demand.var.nonlinear = 'constant';
            m.demand.sigma = 1;
        end
        m.model.endog = false;
        m.demand.alpha = 0.2;
        m.model.beta = [1; -5];
        m.model.x = [5,5];
        m.model.markets = 100;
        m.model.randproducts = false;
        % m.model.products = 5;
        m.model.simulatePrices = false;
        
        m.create();

        m.estDemand.settings.paneltype = paneltype{i};
        results{i} = m.estimate();
    end
    display '*********** Results of LSDV and FE estimation **************'
    disp(results{1})
    disp(results{2})
    display '************************************************************'
    % Loop over cols and rows
%     testtrue(results{1}{'p','Coef'}, results{2}{'p','Coef'}, 1e-2)
%     testtrue(results{1}{'p','Std_err'}, results{2}{'p','Std_err'}, 1e-2)
end


%% Test 8: Nonlinear price, testing findCosts
% Should be improved
display '**********************  Test 8  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.settings.drawmethod = 'quadrature';
m.demand.alpha = 1;
m.demand.sigma = .1;

m.demand.settings.paneltype = 'lsdv';
m.demand.var.nonlinear = 'p';

m.model.endog = true;
m.model.randproducts = true;
m.model.markets = 200;
m.model.products = 10;
m.model.x_vcv = [.2, 1];
m.model.c = 0.9;
m.create();
display(m.model)

results = m.estimate()
testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 2e-1 )

m.findCosts()
testtrue( mean(m.data.c) ,  0.9, 1e-1 )

%% Test 9: Nonlinear price, testing findCosts
% Succeeds at the 1% level
display '**********************  Test 9  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.var.nonlinear = 'p constant';
m.demand.alpha = 1;
m.demand.sigma = [0.1; 1];

m.demand.settings.optimalIV = false;

m.model.endog = false;
m.model.randproducts = false;
m.model.markets = 200;
m.model.products = 10;
m.model.x_vcv = [.2, 1];
m.model.c = 0.9;
m.create();
display(m.model)

results = m.estimate()
testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-2 )