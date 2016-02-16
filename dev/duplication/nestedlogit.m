%% SimMarket mergersim comparison
% Test that Stata and Matlab variable creation are identical
% painkillers is a dataset with mergersim M_* variables included, created
% in Stata by nestedlogit.do.
clear
load painkillers

pk.Ptablets = pk.Ptablets./(pk.cpi/100);
pk.PX = pk.PX./(pk.cpi/100);
    
pk.BL0 = repmat( 2*mean(Estimate.mean(pk.time, pk.X)), size(pk.X)) ;


%% 
% Instruments 

instruments = Estimate.countInstruments(pk, 'time', {'firm', 'form', 'substance'});
instruments(: , [4, 6]) = [];
pk.instruments = instruments + 1;

demand = NLDemand(pk(pk.year < 2009, :));
demand.var.nests = 'form substance';
demand.var.price = 'Ptablets';
demand.var.quantity = 'Xtablets';
demand.var.market = 'time';
demand.var.panel = 'product';
demand.var.marketsize = 'BL0';
demand.var.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];

demand.var.instruments = 'instruments';

%% Demand estimate
demand.settings.estimateMethod = 'gmm';
demand.estimate()
demand.settings.estimateMethod = '2sls';
demand.estimate()
return
%% Demand Settings
disp(demand.settings)
display 'Demand Results:'
disp(demand.results)

%% Elasticities
% Elasticities can be calculated for a market. The market selection has to
% be a single market. It displays the unweighted average elasticities over
% products in the period.

demand.elasticities(pk.year == 2008 & pk.month == 12)
demand.groupElasticities('substance', pk.year == 2008 & pk.month == 12)

%% Market
% Find costs. We set weightedAverages to false to replicate mergersim 1
% behavior.

market = Market(demand);
market.var.firm = 'firm';
market.settings.weightedAverages = false;

market.findCosts(pk.year == 2008 & pk.month == 12)
market.summary()

%% Merger 
% The merger is calculated on a copy of market1, with new ownership.

market2 = copy(market);
market2.var.firm = [];
market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
market2.equilibrium();

display 'Merger results:'
summary(market, market2)

%% Test that results are equal
% All individual prices and costs are within 0.1% of mergersim values

SimMarket.testEqual( pk.Ptablets, pk.M_price, 1e-6);
SimMarket.testEqual( pk.BL0, pk.M_BL0, 1e-3);
SimMarket.testEqual( market.c(pk.year == 2008 & pk.month == 12), ...
    pk.M_costs(pk.year == 2008 & pk.month == 12), 1e-3);
SimMarket.testEqual( market2.p(pk.year == 2008 & pk.month == 12), ...
    pk.M_price2(pk.year == 2008 & pk.month == 12), 1e-4);

