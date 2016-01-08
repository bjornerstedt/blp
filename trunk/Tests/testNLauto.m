%% Test 1: NL auto merger 
% Panel regressions

clear
load('auto.mat')
countrynew = double(T.country);

T.country = double(T.country);
T.segment = double(T.segment);
T.firmname = T.firm;

T.firm = double(T.firm);
T.year2 = T.year .^ 2;

demand = NestedLogitDemand(T);
demand.var.nests = 'segment domestic';
demand.var.quantity = 'qu';
demand.var.marketsize = 'MSIZE';
demand.var.market = 'year country';
demand.var.panel = 'co';
demand.var.exog = ['horsepower fuel width height year year2 '...
    'country2 country3 country4 country5 domestic'];
demand.var.price = 'princ'; % Group shares are automatically endogenous in gls
% endog without as many instruments should give error. No difference in
% ols().
demand.settings.paneltype = 'lsdv';
demand.settings.estimateMethod = 'ols';

demand.init(); % Clean up: Original table, Created vars, demeaned table
demand.estimate();

selection = T.year == 1998 & countrynew== 3 ;

market = Market(demand);
market.var.firm = 'firm';
market.settings.valueShares = false;
market.settings.weightedAverages = false;

market.findCosts(selection);
display Costs:
disp(sum(market.c))

market2 = copy(market);
market2.firm(market2.firm == 15 ) = 26; % Buyer(26) seller(15)

market2.equilibrium(selection);
disp 'Sum new prices:'
disp(sum(market2.p))

result = market.compare(market2);

assert(abs(result{1,'Price2'} - 0.7499176 )<1e-5)
assert(abs(result{1,'PriceCh'} - 0.00265315  )<1e-5)
display '************** Mergersim Test 1 passed ******************'

