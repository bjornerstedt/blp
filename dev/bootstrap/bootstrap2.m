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
mcSim = 50
markets = 30
forecastMarkets = 10

pinc = zeros(3, mcSim, 4);
delta = zeros( 5, mcSim, 4);
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

mergerResult0 = summary(market, market2, 'Selection', m.data.marketid > markets)
pc0 = mergerResult0{:, 'PriceCh'};

%% Standard simulation 
msMarkets = 5;
m.demand.var.instruments = 'c';
m.demand.estimate()
m.findCosts();
market = m.market;

% market.summary()
% market.summary('Selection', m.data.marketid == 1)

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium(m.data.marketid > markets-msMarkets & m.data.marketid <= markets);

mergerResult1 = summary(market, market2, 'Selection', m.data.marketid == markets)
mergerResult1b = summary(market, market2, 'Selection', ...
    m.data.marketid <= markets & m.data.marketid > markets-msMarkets)
pc1 = mergerResult1{:, 'PriceCh'};

%% Expected value simulation 

% predict d
X = mean(m.demand.data.x(m.demand.data.marketid <= markets));
d = m.demand.beta(2).* X + m.demand.beta(3) + ...
    [0;m.demand.results.params.betadummies];

c = accumarray(m.data.productid(m.data.marketid <= markets), ...
    m.market.c(m.data.marketid <= markets), [], @mean );

delta( :, mc, 1) = accumarray(m.data.productid(m.data.marketid > markets),...
    m.data.d(m.data.marketid > markets), [], @mean);
delta( :, mc, 2) = m.data.d(m.data.marketid == markets);
delta( :, mc, 4) = d;

m.demand.d(m.data.marketid == markets+1) = d;
m.market.c(m.data.marketid == markets+1) = c;
market = m.market;

% market.summary()
% market.summary('Selection', m.data.marketid == 1)
market.equilibrium(m.data.marketid == markets+1);

market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium(m.data.marketid == markets+1);

mergerResult2 = summary(market, market2, 'Selection', m.data.marketid == markets+1)
pc2 = mergerResult2{:, 'PriceCh'};

pinc( :, mc, 1) = mergerResult0{:, 'PriceCh'};
pinc( :, mc, 2) = mergerResult1{:, 'PriceCh'};
pinc( :, mc, 3) = mergerResult1b{:, 'PriceCh'};
pinc( :, mc, 4) = mergerResult2{:, 'PriceCh'};


end % mc loop
display 'Mean square forecast errors:'
dmsfe2 = mean(sum((delta(:,:,1) - delta(:,:,2) ) .^2 ))
dmsfe3 = mean(sum((delta(:,:,1) - delta(:,:,3) ) .^2))

msfe2 = mean(sum((pinc(:,:,1) - pinc(:,:,2) ) .^2 ))
msfe3 = mean(sum((pinc(:,:,1) - pinc(:,:,3) ) .^2))
msfe4 = mean(sum((pinc(:,:,1) - pinc(:,:,4) ) .^2))

% True price increase
truePinc = mean(pinc(:,:,1),2)
stderr = std(pinc(:,:,1),0,2)
rmsfe21 = sqrt(mean(((pinc(1,:,1) - pinc(1,:,2) ) .^2 )))
rmsfe31 = sqrt(mean(((pinc(1,:,1) - pinc(1,:,3) ) .^2)))
rmsfe41 = sqrt(mean(((pinc(1,:,1) - pinc(1,:,4) ) .^2)))

% close(hbar);

%% Comparisons

%
% Averages and std deviation
toc
