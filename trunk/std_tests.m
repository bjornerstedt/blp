clear

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2);
testtrue2 = @(x,y)assert(abs((x - y)/y)<2e-2, ...
    'Percentage difference is %f', abs((x - y)/y));
testtrue3 = @(x,y)assert(abs((x - y)/y)<1e-1, ...
    'Percentage difference is %f', abs((x - y)/y));

m = SimMarket();
model = m.model;
model.endog = false;
model.markets = 100;
model.randproducts = false;
model.optimalIV = false;

%% Test 1: NestedLogitDemand
m = SimMarket(model);
m.demand = NestedLogitDemand;
m.init();
sresults = m.simulateDemand()
results = m.estimate()
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'})

testVals = [sresults{1, 'sh'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.016950837696100,5.034299380418873,-1.000887008524414,0.004862493783701];

testSameResults(testVals, correctVals, 1);
%% Test 2: NestedLogitDemand - one nest
m = SimMarket(model);
m.demand = NestedLogitDemand;
m.demand.var.nests = 'type';
m.model.alpha = -2;
m.demand.sigma = 0.5;

m.model.types = 2;

m.init();
m.simulateDemand()
results = m.estimate()
% Test that result is within 0.1% of true value
testtrue1(results{'p','Coef'} , results{'p', 'Beta'})

testVals = [results{'p','Coef'}, results{'lsjg','Coef'}, results{'p','Std_err'}];
correctVals = [-1.996812484535741,0.501859775169641,0.005921296583241];

testSameResults(testVals, correctVals, 2);

%% Test 3: MixedLogitDemand - rc_constant
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'constant';
m.init();
sresults = m.simulateDemand()
results = m.estimate()

testtrue1(results{'p','Coef'} , results{'p', 'Theta'} )

testVals = [sresults{1, 'sh'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.020117449446861,5.052100057072429,-1.002931362162900,0.004913162002386];
testSameResults(testVals, correctVals, 3);

%% Test 4: MixedLogitDemand - rc_x
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'x';
m.init();
sresults = m.simulateDemand()

results = m.estimate()
testVals = [sresults{1, 'sh'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.021868951346931,5.086215427553291,-0.999033679194202,0.005165964169027];
testSameResults(testVals, correctVals, 4);

testtrue(results{'rc_x','Coef'} , m.model.rc_sigma )

%% Test 5: CES MixedLogitDemand
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.settings.ces = true;
%         m.estDemand.settings.robust = 'false';
m.model.endog = false;
m.model.alpha = -4;
m.model.beta = [ 1; 4];
m.model.markets = 200;
m.model.randproducts = false;
m.model.optimalIV = false;
m.demand.var.nonlinear = 'constant';
m.init();
results = m.calculateDemand()
results = m.estimate()
testtrue1(results{'lP','Coef'} , results{'lP', 'Theta'})

%% Test 6: Optimal IV
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'x';

m.model.optimalIV = true;
m.model.endog = false;
m.model.randproducts = true;
m.model.products = 10;
m.init();

results = m.simulateDemand()
results = m.estimate()
testtrue2(results{'p','Coef'} , results{'p', 'Theta'} )

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
        end
        m.model.endog = false;
        m.model.optimalIV = true;
        m.model.alpha = -0.2;
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
    assert(m.approx(results{1}{'p','Coef'}, results{2}{'p','Coef'}))
    assert(m.approx(results{1}{'p','Std_err'}, results{2}{'p','Std_err'}), ...
        'Different Std_err')
end


%% Test 8: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)
model.rc_sigma = .1;

m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.settings.drawmethod = 'quadrature';

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
testtrue3(results{'p','Coef'} , results{'p', 'Theta'} )

m.findCosts(m.estDemand)
testtrue2( mean(m.data.c) ,  0.9 )

% Slightly different demand due to draws
eq2 = m.simulateDemand(m.estDemand)
max( abs((table2array(eq1) - table2array(eq2)) ./table2array(eq1)))

%% Test 9: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)

m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.var.nonlinear = 'p x constant';
m.model.rc_sigma = [0.1; 1; 1];

m.model.optimalIV = false;
m.model.endog = false;
m.model.randproducts = false;
m.model.markets = 100;
m.model.products = 10;
m.model.x_sigma = [.2, 1];
m.model.c = 0.9;
m.init();

eq1 = m.simulateDemand()

results = m.estimate()
testtrue2(results{'p','Coef'} , results{'p', 'Theta'} )