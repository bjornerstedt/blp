clear

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2);
testtrue2 = @(x,y)assert(abs((x - y)/y)<2e-2, ...
    'Percentage difference is %f', abs((x - y)/y));

m = SimMarket();
model = m.model;
model.endog = false;
model.markets = 100;
model.randproducts = false;
model.optimalIV = false;

%% Test: NestedLogitDemand
m = SimMarket(model);
m.demand = NestedLogitDemand;
m.init();
results = m.simulateDemand()
testdiff(results{1, 'sh'} , 0.013209)
testdiff(results{1, 'p'} , 5.0482)
results = m.estimate()
testdiff(results{'p','Coef'} , -0.9996)
testdiff(results{'p','Std_err'} , 0.0049173)
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'})

%% Test: NestedLogitDemand - one nest
m = SimMarket(model);
m.demand = NestedLogitDemand;
m.demand.var.nests = 'type';
m.model.alpha = -2;
m.model.sigma = 0.5;
m.model.types = 2;

m.init();
m.simulateDemand()
results = m.estimate()
testdiff(results{'p','Coef'} , -1.9998)
testdiff(results{'lsjg','Coef'} , 0.49999)
testdiff(results{'p','Std_err'} , 0.0070622)
% Test that result is within 0.1% of true value
testtrue(results{'p','Coef'} , results{'p', 'Beta'})

%% Test: MixedLogitDemand - rc_constant
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.init();
results = m.simulateDemand()
testdiff(results{1, 'sh'} , 0.016533)
testdiff(results{1, 'p'} , 5.0648)
results = m.estimate()
testdiff(results{'p','Coef'} , -1.001)
testdiff(results{'p','Std_err'} , 0.0050055)

testtrue(results{'p','Coef'} , results{'p', 'Theta'} )

%% Test: MixedLogitDemand - rc_x
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.model.nonlinear = 'x';
%m.model.endog = true;
m.init();
results = m.simulateDemand()
testdiff(results{1, 'sh'} , 0.018867)
testdiff(results{1, 'p'} , 5.1041)
results = m.estimate()
testdiff(results{'p','Coef'} , -1.001)
testdiff(results{'p','Std_err'} , 0.0050055)

testtrue(results{'rc_x','Coef'} , m.model.rc_sigma )

%% Test: CES MixedLogitDemand
m = SimMarket(model);
m.demand = MixedLogitDemand;
m.demand.settings.ces = true;
%         m.estDemand.settings.robust = 'false';
m.model.endog = false;
m.model.alpha = -4
m.model.beta = [ 1; 4];
m.model.markets = 200;
m.model.randproducts = false;
m.model.optimalIV = false;
m.init();
results = m.calculateDemand()
results = m.estimate()
testtrue1(results{'lP','Coef'} , results{'lP', 'Theta'})

%% Test: Optimal IV
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

%% Test: LSDV/FE
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


%% Test: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)
model.nonlinear = 'p';
model.rc_sigma = .1;
model.drawmethod = 'quadrature';

m = SimMarket(model);
m.demand = MixedLogitDemand;

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
testtrue2(results{'p','Coef'} , results{'p', 'Theta'} )

m.findCosts(m.estDemand)

% Slightly different demand due to draws
eq2 = m.simulateDemand(m.estDemand)
max( abs((table2array(eq1) - table2array(eq2)) ./table2array(eq1)))

%% Test: Nonlinear price, testing findCosts
% Succeeds at the 2% level (to avoid increasing data size)

m = SimMarket(model);
m.demand = MixedLogitDemand;
m.model.nonlinear = 'p x constant';
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