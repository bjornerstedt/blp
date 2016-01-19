% Test that Stata and Matlab variable creation are identical

clear
pkS = readtable('test.csv','TreatAsEmpty','NA');
load painkillers

pk.Ptablets = pk.Ptablets./(pk.cpi/100);
pk.PX = pk.PX./(pk.cpi/100);
    
[~,~, index] = unique(pk.time);    
pk.BL0 = repmat( 2*mean(accumarray(index, pk.X)), size(pk.X)) ;

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

demand.estimate()
demand.settings
demand.results

market = Market(demand);
market.var.firm = 'firm';
% mergersim 1 has unweighted averages
market.settings.weightedAverages = false;

market.findCosts(pk.year == 2008 & pk.month == 12)

market.results.findCosts
market.summary()

market2 = copy(market);
market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
market2.equilibrium();
compare(market, market2)

%%%% Test that results are equal
SimMarket.testEqual( pk.Ptablets, pk.M_price, 1e-6);
SimMarket.testEqual( pk.BL0, pk.M_BL0, 1e-3);
SimMarket.testEqual( market.c(pk.year == 2008 & pk.month == 12), ...
    pk.M_costs(pk.year == 2008 & pk.month == 12), 1e-3);
SimMarket.testEqual( market2.p(pk.year == 2008 & pk.month == 12), ...
    pk.M_price2(pk.year == 2008 & pk.month == 12), 1e-4);
%%%%
