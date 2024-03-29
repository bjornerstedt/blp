clear
% Simulated Testing

%% Test 1: NLDemand
display '**********************  Test 1  *************************'
m = SimMarket();
m.model.markets = 200;
m.demand = NLDemand;
m.demand.alpha = 1;

sresults = m.create();
display(m.model)
results = m.demand.estimate()
% Test that result is within 0.1% of true value
SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2)

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.011584066958451,5.031660735481728,-1.002459854394107,0.003067662778767];

SimMarket.testSame(testVals, correctVals, 1);

m.findCosts()
SimMarket.testEqual( mean(m.data.c) ,  m.model.c, 1e-2 )

%% Test 2: NLDemand - one nest CES demand
display '**********************  Test 2  *************************'
m = SimMarket();
m.model.markets = 200;
m.demand = NLDemand;
m.demand.settings.ces = true;
m.demand.var.nests = 'type';
m.demand.alpha = 2;
m.demand.sigma = 0.5;
m.model.typeList = [1,2,1,2,1];
m.create();
display(m.model)

results = m.demand.estimate()
% Test that result is within 0.1% of true value
SimMarket.testEqual(results{'lP','Coef'} , results{'lP', 'True_val'}, 1e-2)

testVals = [results{'lP','Coef'}, results{'lsjg','Coef'}, results{'lP','Std_err'}];
correctVals = [-1.997982777266492,0.506308608342166,0.014891089529195];

SimMarket.testSame(testVals, correctVals, 2);
m.findCosts()
SimMarket.testEqual( mean(m.data.c) ,  m.model.c, 1e-2 )

%% Test 3: NLDemand - two level nests
display '**********************  Test 3 *************************'
m = SimMarket();
m.demand = NLDemand;
m.demand.var.nests = 'type1 type2';
m.demand.alpha = 2;
m.demand.sigma = [0.6, 0.3];
m.model.types = [2, 3];
m.model.products = 9;
m.create();

display(m.model)

results = m.demand.estimate()
% Test that result is within 0.1% of true value
SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2)

testVals = [results{'p','Coef'}, results{'lsjh','Coef'}, results{'p','Std_err'}];
correctVals = [-1.994459335261852,0.601735881106264,0.004873164058653];

SimMarket.testSame(testVals, correctVals, 3);

%% Test 4: RCDemand - rc_constant
% Estimate without instruments, with exogeneous p.
display '**********************  Test 4  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.var.nonlinear = 'constant';
sresults = m.create();
display(m.model)

results = m.demand.estimate()
% Check that estimate different than initial guess:
assert(abs(m.demand.results.sigma0 - m.demand.sigma) > 0.3)

SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2 )
SimMarket.testEqual(results{'rc_constant','Coef'} , results{'rc_constant', 'True_val'}, 1e-1 )

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.027209827330165,4.980964755245374,-1.002396191815819,0.004427962469071];
SimMarket.testSame(testVals, correctVals, 4);

%% Test 5: RCDemand - rc_x
display '**********************  Test 5  *************************'
% demchoice = {RCDemand, RCDemand2};
for i = 0:0
    m = SimMarket();
    m.demand = RCDemand;
    m.demand.config.compiled = logical(i);
    m.demand.alpha = 1;
    m.demand.sigma = 1;
    m.demand.var.nonlinear = 'x';
    sresults = m.create();
    display(m.model)
    
    results = m.demand.estimate()
    testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
    correctVals = [0.030538586199379,4.980964755245374,-1.005205788003259,0.004425086528973];
    SimMarket.testSame(testVals, correctVals, 5);
    
    SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2 )
    SimMarket.testEqual(results{'rc_x','Coef'} , results{'rc_x', 'True_val'}, 1e-1 )
end
%% Test 6: CES RCDemand, exogenous
display '**********************  Test 6  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true;
m.demand.alpha = 3;
m.demand.sigma = .5;
m.model.beta = [ 1, 4];
m.model.markets = 200;
m.demand.var.nonlinear = 'p';
% m.model.endog = true;
% m.model.randomProducts = false;
% m.model.pricesFromCosts = false;

m.create();
display(m.model)

results = m.demand.estimate()
testVals = [results{'lP', 'Coef'}, results{'rc_lP', 'Coef'}];
correctVals = [-3.012747977629290,0.511180973525220];
SimMarket.testSame(testVals, correctVals, 6);

SimMarket.testEqual(results{'lP', 'Coef'} , results{'lP', 'True_val'}, 1e-2)
SimMarket.testEqual(results{'rc_lP', 'Coef'} , results{'rc_lP', 'True_val'}, 5e-2)

%% Test 7: Optimal IV
display '**********************  Test 7  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.var.nonlinear = 'x';
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.settings.optimalIV = true;
% m.demand.settings.ces = true; % This works

m.model.endog = true;
m.model.randomProducts = true;
m.model.products = 10;
m.create();
display(m.model)
% Ad-hoc to keep quadratic instruments
m.demand.data.nprod2 = m.data.nprod.^2;
m.demand.var.instruments = 'nprod nprod2'
results = m.demand.estimate()

testVals = [results{'p','Coef'}, results{'rc_x','Coef'}];
correctVals = [-0.970459257509805,1.045491435731395];

SimMarket.testSame(testVals, correctVals, 7);

SimMarket.testEqual(results{'rc_x','Coef'} , results{'rc_x', 'True_val'}, 5e-2 )

%% Test 8: LSDV/FE
% Four tests, comparing LSDV/FE for NL/FE
display '**********************  Test 8  *************************'
results = cell(2,1);
for d = 1:2
    paneltype = {'lsdv', 'fe'};
    for i = 1:2
        m = SimMarket();
       if d == 1
            m.demand = NLDemand;
        else
            m.demand = RCDemand;
            m.demand.var.nonlinear = 'p';
            m.demand.sigma = .2;
        end
        m.demand.alpha = 1;
%         m.demand.var.instruments = 'c nprod nprod2';
%         m.model.beta = [1; -5];
%         m.model.x = [1,0];
%         m.model.x_vcv = [.2,1];
%         m.model.markets = 100;
        m.model.endog = true;
        m.model.randomProducts = true;
%        m.model.pricesFromCosts = false;
        m.demand.settings.paneltype = paneltype{i};
        
        m.create();
% Ad-hoc to keep quadratic instruments
m.demand.data.nprod2 = m.data.nprod.^2;
m.demand.var.instruments = 'c nprod nprod2'

        results{i} = m.demand.estimate();
        display(results{1})
        display(results{2})
        market = Market(m.demand);
        market.var.firm = 'productid';
       
        display 'Test that mean calculated costs are close to actual'
        market.findCosts( );
        SimMarket.testEqual(mean(m.data.c), mean(market.c), 1e-1)
        
        display 'Test that mean equilibrium prices are close to starting vals'
        market.equilibrium( )
        SimMarket.testEqual(mean(m.data.p), mean(market.p), 1e-4)
    end
    display '*********** Results of LSDV and FE estimation **************'
    display '************************************************************'
    % Loop over cols and rows
    SimMarket.testEqual(results{1}{'p','Coef'}, results{2}{'p','Coef'}, 1e-4)
    SimMarket.testEqual(results{1}{'p','Std_err'}, results{2}{'p','Std_err'}, 1e-4)
end

if false
%% Test 9: Exogenous nonlinear price constant, testing findCosts
display '**********************  Test 9  *************************'
for tst = 0:1
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true; 
m.demand.config.compiled = false;
% m.demand.config.tolerance = 1e-12;
m.demand.settings.drawMethod = 'quadrature';
m.demand.config.guessDelta = false;
m.demand.config.fpmaxit = 3000; 
m.demand.var.exog = 'x1 x2';
m.demand.alpha = 1;
if logical(tst)
    m.demand.sigma = [1; .2];
    m.demand.var.nonlinear = {'x1 x2', 'normal'};
else
    m.demand.sigma = [.2; 1];
    m.demand.var.nonlinear = {'x2 x1', 'normal'};
end
m.demand.settings.nind = 1000;
%  m.demand.settings.marketDraws = true;
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

m.findCosts()
meanCosts = mean(m.data.c)
% SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )
end
end 
%% Test 10: Nonlinear price, testing findCosts
display '**********************  Test 10  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.settings.drawMethod = 'quadrature';
m.demand.alpha = 1;
m.demand.sigma = [0.2; 1];

% m.demand.settings.ces = true;
m.demand.settings.paneltype = 'lsdv';
m.demand.var.nonlinear = 'p x';
m.demand.settings.optimalIV = true;
m.model.gamma = 1;
m.model.endog = true;
m.model.randomProducts = true;
m.model.markets = 200;
m.demand.config.guessDelta = false;
m.model.products = 5;
m.create();
display(m.model)
% Ad-hoc to keep quadratic instruments
m.demand.data.nprod2 = m.data.nprod.^2;
m.demand.var.instruments = 'nprod nprod2 c'

results = m.demand.estimate()
SimMarket.testEqual(results{'rc_p','Coef'} , results{'rc_p', 'True_val'}, 2e-1 )

m.findCosts()
meanCosts = mean(m.data.c)
SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )

%% Test 11: Nonlinear price with CES, multiple distributions
display '**********************  Test 11  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.settings.drawMethod = 'quadrature';
m.demand.alpha = 1;
m.demand.sigma = [0.2; .1];

m.demand.var.nonlinear = {{'p', 'uniform'}, {'x'}};

m.demand.settings.paneltype = 'lsdv';
% m.demand.var.instruments = 'nprod nprod2 c';
m.demand.settings.optimalIV = true;
m.model.gamma = 1;
m.model.endog = true;
m.model.randomProducts = true;
m.model.markets = 200;
m.demand.config.guessDelta = false;
m.model.products = 5;
m.create();
display(m.model)
% Ad-hoc to keep quadratic instruments
m.demand.data.nprod2 = m.data.nprod.^2;
m.demand.var.instruments = 'nprod nprod2 c'


results = m.demand.estimate()
SimMarket.testEqual(results{'rc_p','Coef'} , results{'rc_p', 'True_val'}, 5e-2 )
sigma0 = m.demand.results.sigma0 

