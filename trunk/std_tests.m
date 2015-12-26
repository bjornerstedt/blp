clear

% Put these last to avoid testi
testSameResults = @(x,y,i)true;
testtrue = @(x,y,z)true;

testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);
testtrue = @(x,y,z)assert(abs((x - y)/y)<z, 'Percentage diff is %f', abs((x - y)/y));

%% Test 1: NestedLogitDemand
display '**********************  Test 1  *************************'
m = SimMarket();
m.demand = NestedLogitDemand;
m.demand.alpha = 1;
sresults = m.create();
display(m.model)
results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'}, 1e-2)

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.012851335425658,5.201998556464783,-1.002853220422230,0.004537307509000];

testSameResults(testVals, correctVals, 1);
%% Test 2: NestedLogitDemand - one nest
display '**********************  Test 2  *************************'
m = SimMarket();
m.demand = NestedLogitDemand;
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
correctVals = [-2.003006191822202,0.499958011734008,0.006830329622917];

testSameResults(testVals, correctVals, 2);

%% Test 3: MixedLogitDemand - rc_constant
% Estimate without instruments, with exogeneous p.
display '**********************  Test 3  *************************'
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.alpha = 1;
m.demand.rc_sigma = 1;
m.demand.var.nonlinear = 'constant';
sresults = m.create();
display(m.model)

results = m.estimate()
% Check that estimate different than initial guess:
assert(abs(m.estDemand.results.rc_sigma0 - m.estDemand.rc_sigma) > 0.5)

testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-2 )

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.016785830391304,5.223475579061942,-1.001934217196106,0.004867296895755];
testSameResults(testVals, correctVals, 3);

%% Test 4: MixedLogitDemand - rc_x
display '**********************  Test 4  *************************'
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.alpha = 1;
m.demand.rc_sigma = 1;
m.demand.var.nonlinear = 'x';
sresults = m.create();
display(m.model)

results = m.estimate()
testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.017205956573730,5.252125999940534,-1.003194332212289,0.004661671210179];
testSameResults(testVals, correctVals, 4);

%testtrue(results{'rc_x','Coef'} , m.demand.rc_sigma, 1e-1 )

%% Test 5: CES MixedLogitDemand
display '**********************  Test 5  *************************'
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.settings.ces = true;
%         m.estDemand.settings.robust = 'false';
m.model.endog = false;
m.demand.alpha = 4;
m.demand.rc_sigma = 1;
m.model.beta = [ 1; 4];
m.model.markets = 200;
m.model.randproducts = false;
m.demand.var.nonlinear = 'constant';
m.model.simulatePrices = false;

m.create();
display(m.model)

results = m.estimate()
testVals = [results{'lP', 'Coef'}, results{'rc_constant', 'Coef'}];
correctVals = [-3.993081383179611,1.000792980033876];
testSameResults(testVals, correctVals, 5);

testtrue(results{'lP', 'Coef'} , results{'lP', 'Theta'}, 2e-2)

%% Test 6: Optimal IV
display '**********************  Test 6  *************************'
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'x';
m.demand.alpha = 1;
m.demand.rc_sigma = 1;
m.demand.settings.optimalIV = true;

m.model.endog = true;
m.model.randproducts = true;
m.model.products = 10;
m.create();
display(m.model)

results = m.estimate()

testVals = [results{'p','Coef'}, results{'rc_x','Coef'}];
correctVals = [-0.661577893863166,0.961761518545162];
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
            m.demand = NestedLogitDemand;
        else
            m.demand = MixedLogitDemand;
            m.demand.var.nonlinear = 'constant';
            m.demand.rc_sigma = 1;
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
% Succeeds at the 2% level (to avoid increasing data size)
display '**********************  Test 8  *************************'

m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.settings.drawmethod = 'quadrature';
m.demand.alpha = 1;
m.demand.rc_sigma = .1;

m.demand.settings.paneltype = 'lsdv';
m.demand.var.nonlinear = 'p';

m.model.endog = true;
m.model.randproducts = true;
m.model.markets = 100;
m.model.products = 10;
m.model.x_sigma = [.2, 1];
m.model.c = 0.9;
m.create();
display(m.model)

results = m.estimate()
testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-1 )

m.findCosts(m.estDemand)
testtrue( mean(m.data.c) ,  0.9, 1e-1 )

%% Test 9: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)
display '**********************  Test 9  *************************'

m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'p constant';
m.demand.alpha = 1;
m.demand.rc_sigma = [0.1; 1];

m.demand.settings.optimalIV = false;

m.model.endog = false;
m.model.randproducts = false;
m.model.markets = 100;
m.model.products = 10;
m.model.x_sigma = [.2, 1];
m.model.c = 0.9;
m.create();
display(m.model)

results = m.estimate()
testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 2e-2 )