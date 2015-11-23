% Relevant tests for modification of MixedLogitDemandNew
clear

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2);

m = SimMarket();
model = m.model;
model.endog = false;
model.markets = 100;
model.randproducts = false;
model.optimalIV = false;

%% Test: MixedLogitDemand
m0 = SimMarket(model);
m0.demand = MixedLogitDemand;
m0.init();
m0.calculateDemand();
d0 = m0.estDemand;
d0.init();
expmu = d0.nlpart(m0.model.rc_sigma, d0.vx);
edelta = d0.share.s ./ d0.share.s0;
edelta = d0.finddelta(expmu, edelta, 1e-12);
% edelta(m0.data.marketid == 1)
[~,indsh] = d0.sharecalc( edelta, expmu);
% indsh(m0.data.marketid == 1, 1:8)

sj = d0.deltaJacobian(m0.model.rc_sigma, edelta);
sj(m0.data.marketid == 1)

%% Test: MixedLogitDemandNew
m = SimMarket(model);
m.demand = MixedLogitDemandNew;
m.init();
m.calculateDemand();
d = m.estDemand;
d.init();
d.nlpart(m.model.rc_sigma, d.vx);
dp = d.period{1};
% dp.expmu(:,1:8)

edelta = d.share.s ./ d.share.s0;
edelta = edelta(m.data.marketid == 1);
edelta = dp.finddelta(m.model.rc_sigma, edelta, 1e-12);
[~,indsh] = dp.sharecalc( edelta);
% indsh(:, 1:8)
sj = dp.deltaJacobian(m.model.rc_sigma, edelta)

return




results = m.simulateDemand()
testdiff(results{1, 'sh'} , 0.016533)
testdiff(results{1, 'p'} , 5.0648)
results = m.estimate()
testdiff(results{'p','Coef'} , -1.001)
testdiff(results{'p','Std_err'} , 0.0050055)

testtrue(results{'p','Coef'} , m.model.beta(1) )

%% Test: CES MixedLogitDemandNew
m = SimMarket(model);
m.demand = MixedLogitDemandNew;
m.model.ces = true;
%         m.estDemand.settings.robust = 'false';
m.model.endog = false;
m.model.beta = [-4; 1; 4];
m.model.markets = 200;
m.model.randproducts = false;
m.model.optimalIV = false;
m.init();
results = m.calculateDemand()
result = m.estimate()
testtrue1(result{'lP','Coef'} , m.model.beta(1))

%% Test: Optimal IV
m = SimMarket(model);
m.demand = MixedLogitDemandNew;
m.demand.var.nonlinear = 'x';

m.model.optimalIV = true;
m.model.endog = false;
m.model.randproducts = true;
m.model.products = 10;
m.init();

results = m.simulateDemand()
results = m.estimate()
testtrue1(results{'p','Coef'} , m.model.beta(1) )

%% Test: LSDV/FE
% Four tests, comparing LSDV/FE for NL/FE
results = cell(2,1);
paneltype = {'lsdv', 'fe'};
for i = 1:2
    m = SimMarket();
    m.demand = MixedLogitDemandNew;
    m.model.endog = false;
    m.model.optimalIV = true;
    m.model.beta = [-.2; 1; -5];
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

