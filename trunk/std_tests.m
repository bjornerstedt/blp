clear

% Put these last to avoid testi
testSameResults = @(x,y,i)true;
testtrue = @(x,y,z)true;


testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);
testtrue = @(x,y,z)assert(abs((x - y)/y)<z, 'Percentage diff is %f', abs((x - y)/y));

%% Test 1: NestedLogitDemand
m = SimMarket();
m.demand = NestedLogitDemand;
m.demand.alpha = 1;
m.init();
sresults = m.simulateDemand()
results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'}, 1e-3)

testVals = [sresults{1, 'sh'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.016950837696100,5.034299380418873,-1.000887008524414,0.004862493783701];

testSameResults(testVals, correctVals, 1);
%% Test 2: NestedLogitDemand - one nest
m = SimMarket();
m.demand = NestedLogitDemand;
m.demand.var.nests = 'type';
m.demand.alpha = 2;
m.demand.sigma = 0.5;

m.model.types = 2;

m.init();
m.simulateDemand()
results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'}, 1e-2)

testVals = [results{'p','Coef'}, results{'lsjg','Coef'}, results{'p','Std_err'}];
correctVals = [-1.996812484535741,0.501859775169641,0.005921296583241];

testSameResults(testVals, correctVals, 2);

%% Test 3: MixedLogitDemand - rc_constant
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.alpha = 1;
m.demand.rc_sigma = 1;
m.demand.var.nonlinear = 'constant';
m.init();
sresults = m.simulateDemand()
results = m.estimate()

testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-2 )

testVals = [sresults{1, 'sh'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.020117449446861,5.052100057072429,-1.008306669278507,0.005101784438227];
testSameResults(testVals, correctVals, 3);

%% Test 4: MixedLogitDemand - rc_x
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.alpha = 1;
m.demand.rc_sigma = 1;
m.demand.var.nonlinear = 'x';
m.init();
sresults = m.simulateDemand()

results = m.estimate()
testVals = [sresults{1, 'sh'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.021868951346931,5.086215427553291,-1.008741759875677,0.005362773381715];
testSameResults(testVals, correctVals, 4);

%testtrue(results{'rc_x','Coef'} , m.demand.rc_sigma, 1e-1 )

%% Test 5: CES MixedLogitDemand
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
m.init();
results = m.calculateDemand()
results = m.estimate()

testVals = [results{'lP','Coef'}, results{'rc_constant','Coef'}];
correctVals = [-3.929836828254647,0.760600817538613];
testSameResults(testVals, correctVals, 5);

testtrue(results{'lP','Coef'} , results{'lP', 'Theta'}, 2e-2)

%% Test 6: Optimal IV
m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'x';
m.demand.alpha = 1;
m.demand.rc_sigma = 1;

m.model.optimalIV = true;
m.model.endog = false;
m.model.randproducts = true;
m.model.products = 10;
m.init();
%m.demand.settings.optimalIV = true;

results = m.simulateDemand()
results = m.estimate()

testVals = [results{'p','Coef'}, results{'rc_x','Coef'}];
correctVals = [-0.926700943042958,0.251569943194395];
testSameResults(testVals, correctVals, 6);

testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-1 )

%% Test 7: LSDV/FE
% Four tests, comparing LSDV/FE for NL/FE
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
        m.model.optimalIV = true;
        m.demand.alpha = 0.2;
        m.model.beta = [1; -5];
        m.model.x = [5,5];
        m.model.markets = 100;
        m.model.randproducts = false;
        % m.model.products = 5;
        
        m.init();
        m.calculateDemand()
        %m.simulateDemand()
        m.estDemand.settings.paneltype = paneltype{i};
        results{i} = m.estimate();
    end
    display '*********** Results of LSDV and FE estimation **************'
    disp(results{1})
    disp(results{2})
    display '************************************************************'
    % Loop over cols and rows
    testtrue(results{1}{'p','Coef'}, results{2}{'p','Coef'}, 1e-2)
    testtrue(results{1}{'p','Std_err'}, results{2}{'p','Std_err'}, 1e-2)
end


%% Test 8: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)

m = SimMarket();
m.demand = MixedLogitDemand;
m.demand.settings.drawmethod = 'quadrature';
m.demand.alpha = 1;
m.demand.rc_sigma = .1;

m.demand.settings.paneltype = 'lsdv';
m.demand.var.nonlinear = 'p';
m.model.optimalIV = false;
m.model.endog = true;
m.model.randproducts = true;
m.model.markets = 100;
m.model.products = 10;
m.model.x_sigma = [.2, 1];
m.model.c = 0.9;
m.init();

eq1 = m.simulateDemand()

results = m.estimate()
testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 1e-1 )

m.findCosts(m.estDemand)
testtrue( mean(m.data.c) ,  0.9, 1e-1 )

% Slightly different demand due to draws
eq2 = m.simulateDemand(m.estDemand)
max( abs((table2array(eq1) - table2array(eq2)) ./table2array(eq1)))

%% Test 9: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)

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
m.init();

eq1 = m.simulateDemand()

results = m.estimate()
testtrue(results{'p','Coef'} , results{'p', 'Theta'}, 2e-2 )