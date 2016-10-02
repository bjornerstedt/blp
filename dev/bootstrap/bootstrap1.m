%% Merger with estimated demand
if true
display '**********************   Generate sample  *************************'
m = SimMarket();
m.model.markets = 50;
m.demand = NLDemand;
m.demand.alpha = .5;
m.model.endog = true;
m.model.randomProducts = true;
m.model.firm = [1,1,2,2,3];
m.model.beta(2) = 1;

% [mk, beta0] = m.calibrate(.5)
m.create()
sum(mk{:,2})

m.demand.data.firm2 = m.data.firm;
m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;

display(m.model)
results = m.demand.estimate()


m.findCosts();
market = m.market;

market.summary()
market.summary('Selection', m.data.marketid == 1)


market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium();

mergerResult1 = summary(market, market2);
display 'Price change: '
pc = mergerResult1{1, 'PriceCh'}
end

display '**********************  Bootstraps  *************************'

for r = 1:1

%% Parametric bootstraps alpha / alpha & beta0
% 
% Here we estimate and use repeated draws of alpha to calculate expected
% price increase. Use one loop to do both parametric and nonparametric
% bootstraps.

% Estimate
% Loop taking random draws of
% alpha
% alpha and beta0

%% Non-parametric bootstrap 

% Take bootstrap draw of data
y = datasample(m.data, size(m.data, 1)) 

% estimate

end

%% Comparisons
%
% Averages and std deviation