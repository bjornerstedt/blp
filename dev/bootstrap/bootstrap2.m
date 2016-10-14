%% Merger with estimated demand
% This script compares the simulated effects of a merger compared with the
% actual difference. 

% In the forecast period we study the simulated effect using:
% 1) the true parameters 
% 2) a parametric bootstrap of alpha in the last period. 
% 3) a parametric bootstrap of expected values. 
% 4) a nonparametric bootstrap estimate and simulation

tic
% SimMarket.randDraws(22);
mcSim = 2
markets = 25
forecastMarkets = 1

pinc = zeros(mcSim, 4);
% hbar = parfor_progressbar(mcSim,'Computing...');
for mc = 1:mcSim
display(sprintf('*******************   Simulation %d of %d  **********************', mc, mcSim));

m = SimMarket();
m.model.markets = markets + forecastMarkets;
m.demand = NLDemand;
m.demand.alpha = .5;
m.model.endog = true;
% m.model.randomProducts = true;
m.model.firm = [1,1,2,2,3];
m.model.beta(2) = 1;
m.config.initRandomDraws = false;

% [mk, beta0] = m.calibrate(.5)
mk = m.create();
% sum(mk{:,2})

% Adding a field to data means both demand and market must be updated
m.demand.data.firm2 = m.data.firm;
m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;
m.market.demand = m.demand;

% display(m.model)
%% True outcome
market = m.market;

% market.summary()
% market.summary('Selection', m.data.marketid == 1)

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium(m.data.marketid > markets);

mergerResult0 = summary(market, market2)
pc0 = mergerResult0{1, 'PriceCh'};

%% Standard simulation 
m.demand.var.instruments = 'c';
m.demand.estimate()
m.findCosts();
market = m.market;

% market.summary()
% market.summary('Selection', m.data.marketid == 1)

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium(m.data.marketid == markets);

mergerResult1 = summary(market, market2)

%% Expected value simulation 

% predict d
X = mean(m.demand.data.x(m.demand.data.marketid <= markets));
d = m.demand.beta(2).* X + m.demand.beta(3) + ...
    [0;m.demand.results.params.betadummies];

d
m.data.d(m.data.marketid > markets)
m.data.d(m.data.marketid == markets)





return
m.findCosts();
market = m.market;

% market.summary()
% market.summary('Selection', m.data.marketid == 1)

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium(m.data.marketid == markets);

mergerResult1 = summary(market, market2)

return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vcv = m.demand.results.params.varcovar;

% High correlation:
% corrbeta = vcv{'p', 'constant'} / sqrt(vcv{'p', 'p'}) / sqrt(vcv{'constant', 'constant'});

%% Parametric bootstraps alpha / alpha & beta0
% 
% Here we estimate and use repeated draws of alpha to calculate expected
% price increase. Use one loop to do both parametric and nonparametric
% bootstraps.
%
% Estimate
% Loop taking random draws of
% alpha 
% alpha x and beta0

sel = {[1], [1,2,3]};

bootpars = parametricBootstrap(m.demand.results.estimate{sel{1},'Coef'}, vcv{sel{1},sel{1}}, bootreps) ;
data = m.demand.data;
demand = copy(m.demand);
pc = zeros(bootreps, 3);
for r = 1:bootreps
if true
    
demand.beta(1) = bootpars(r, :)';
demand.alpha = -demand.beta(1);
demand.calibrate();

market = Market(demand);
market.config.quietly = true;

market.var.firm = 'firm';
market.findCosts();

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium();

mergerResult1 = summary(market, market2);
pc(r,1) = demand.alpha;
pc(r, 2) = mergerResult1{1, 'PriceCh'};
delete(market2)

end

if false
%% Non-parametric bootstrap 
newdemand = copy(m.demand);
newdemand.data = bootstrap(data, 'marketid');
newdemand.estimate();

market = Market(newdemand);

market.var.firm = 'firm';
market.findCosts();

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium();

mergerResult1 = summary(market, market2);
pc(r, 3) = mergerResult1{1, 'PriceCh'};
delete(market2)
end


end % pb loop
pinc(mc,:) = [pc0, mean(pc,1)];
hbar.iterate(1);
end % mc loop
close(hbar);

%% Comparisons
pinc
%
% Averages and std deviation
display Averages:
mean(pinc)
toc
