clear
% Simulated Testing

testSameResults = @(x,y,i)assert(all( abs(x - y) < 10e-6), 'Test %d failed' , i);
% Uncomment to avoid testing
% testSameResults = @(x,y,i)true;


%% Test 1: NLDemand
display '**********************  Test 1  *************************'
m = SimMarket();
m.demand = NLDemand;
m.demand.alpha = 1;

sresults = m.create();
display(m.model)
results = m.estimate()
% Test that result is within 0.1% of true value
SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2)

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.021574061753928,4.980964755245374,-1.003422555598600,0.004304018905541];

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
SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2)

testVals = [results{'p','Coef'}, results{'lsjg','Coef'}, results{'p','Std_err'}];
correctVals = [-2.003706029185489,0.499863507432290,0.005375770177673];

testSameResults(testVals, correctVals, 2);

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

results = m.estimate()
% Test that result is within 0.1% of true value
SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2)

testVals = [results{'p','Coef'}, results{'lsjh','Coef'}, results{'p','Std_err'}];
correctVals = [-1.994459335261852,0.601735881106264,0.004873164058653];

testSameResults(testVals, correctVals, 3);

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

results = m.estimate()
% Check that estimate different than initial guess:
assert(abs(m.estDemand.results.sigma0 - m.estDemand.sigma) > 0.3)

SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2 )

testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.027209827330165,4.980964755245374,-1.002058277997551,0.004437443811257];
testSameResults(testVals, correctVals, 4);

%% Test 5: RCDemand - rc_x
display '**********************  Test 5  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.alpha = 1;
m.demand.sigma = 1;
m.demand.var.nonlinear = 'x';
sresults = m.create();
display(m.model)

results = m.estimate()
testVals = [sresults{1, 'q'}, sresults{1, 'p'}, results{'p','Coef'}, results{'p','Std_err'}];
correctVals = [0.030538586199379,4.980964755245374,-1.004698788631557,0.004361554728947];
testSameResults(testVals, correctVals, 5);

SimMarket.testEqual(results{'p','Coef'} , results{'p', 'True_val'}, 1e-2 )
SimMarket.testEqual(results{'rc_x','Coef'} , results{'rc_x', 'True_val'}, 2e-2 )

%% Test 6: CES RCDemand
display '**********************  Test 6  *************************'
m = SimMarket();
m.demand = RCDemand;
m.demand.settings.ces = true;
%         m.estDemand.settings.robust = 'false';
m.demand.alpha = 4;
m.demand.sigma = 1;
m.model.beta = [ 1; 4];
% m.model.markets = 200;
m.demand.var.nonlinear = 'constant';
% m.model.endog = true;
% m.model.randomProducts = false;
% m.model.pricesFromCosts = false;

m.create();
display(m.model)

results = m.estimate()
testVals = [results{'lP', 'Coef'}, results{'rc_constant', 'Coef'}];
correctVals = [-4.012600987612641,0.973319662065424];
testSameResults(testVals, correctVals, 6);

SimMarket.testEqual(results{'lP', 'Coef'} , results{'lP', 'True_val'}, 1e-2)

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
        m.demand.var.instruments = 'c nprod nprod2';
%         m.model.beta = [1; -5];
%         m.model.x = [1,0];
%         m.model.x_vcv = [.2,1];
%         m.model.markets = 100;
        m.model.endog = true;
        m.model.randomProducts = true;
%        m.model.pricesFromCosts = false;
        
        m.create();

        m.estDemand.settings.paneltype = paneltype{i};
        results{i} = m.estimate();
        display(results{1})
        display(results{2})
        market = Market(m.estDemand);
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

%% Test 9: Exogenous nonlinear price constant, testing findCosts
display '**********************  Test 9  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.var.nonlinear = 'p x';
m.demand.alpha = 1;
m.demand.sigma = [0.5; 1];
m.demand.settings.nind = 500;

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

%% Test 10: Nonlinear price, testing findCosts
display '**********************  Test 10  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.settings.drawmethod = 'quadrature';
m.demand.alpha = 1;
m.demand.sigma = [0.2; 1];

% m.demand.settings.ces = true;
m.demand.settings.paneltype = 'lsdv';
m.demand.var.nonlinear = 'p x';
m.demand.var.instruments = 'nprod nprod2 c';
m.demand.settings.optimalIV = true;
m.model.gamma = 1;
m.model.endog = true;
m.model.randomProducts = true;
m.model.markets = 200;
m.demand.config.guessdelta = false;
m.model.products = 5;
m.create();
display(m.model)

results = m.estimate()
SimMarket.testEqual(results{'rc_p','Coef'} , results{'rc_p', 'True_val'}, 2e-1 )

m.findCosts()
meanCosts = mean(m.data.c)
SimMarket.testEqual( meanCosts ,  m.model.c, 1e-2 )

%% Test 11: Nonlinear price with CES, multiple distributions
display '**********************  Test 11  *************************'

m = SimMarket();
m.demand = RCDemand;
m.demand.settings.drawmethod = 'quadrature';
m.demand.alpha = 1;
m.demand.sigma = [0.2; 1];

m.demand.var.nonlinear = {{'p', 'uniform'}, {'x'}};

m.demand.settings.paneltype = 'lsdv';
m.demand.var.instruments = 'nprod nprod2 c';
m.demand.settings.optimalIV = true;
m.model.gamma = 1;
m.model.endog = true;
m.model.randomProducts = true;
m.model.markets = 200;
m.demand.config.guessdelta = false;
m.model.products = 5;
m.create();
display(m.model)

results = m.estimate()
SimMarket.testEqual(results{'rc_p','Coef'} , results{'rc_p', 'True_val'}, 5e-2 )

