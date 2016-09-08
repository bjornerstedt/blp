%% Unbalanced nests
% Does

demand = NLDemand();
demand.alpha = 0.5;
demand.sigma = 0.5;
demand.var.nests = 'type';

m3 = SimMarket();
m3.demand = demand;
m3.model.typeList = [1,1,1,1,2];
%m3.model.typeList = [1,2,1,2,1];
m3.model.typeList = [1,1,3,1,2];

m3.market = Market;

m3.model.firm = [1,1,1,2,2];
% m3.model.firm = [1,2,3,4,5];

m3.model.endog = true;
m3.model.randomProducts = true;

m3.model.gamma = 1;
m3.create();
dt3 = m3.data;

%%
% As before, we can estimate the model by first creating a demand object and
% then specifying a set of parameters. Here we wait with associating the
% dataset dt3 with the demand, as we want to add a set of instruments to
% the dataset first.
demand = NLDemand;

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketSize = 'constant';
demand.var.exog = 'x';

%%
% We can create an array of count instruments by using the utility method
% |Estimate.countInstruments()|. It creates count instruments by market as
% specified in the 'marketid' column, here creating counts by market and
% all combinations of firms and type. Note that a column of a Matlab table
% can be an array instead of a column vector. Here the count instruments
% are all put in the column |dt3.inst| of the
dt3.inst = Estimate.countInstruments(dt3, 'marketid', {'firm', 'type'});
demand.var.instruments = 'inst';
demand.data = dt3;
%%
% To estimate this model, the nesting variable |type| has to be specified.
% The same count instruments as above are used. Note that the 2SLS FE panel
% estimate will have price and log group shares as endogenous variables.

demand.var.nests = 'type';
demand.settings.paneltype = 'lsdv';
result = demand.estimate()

%% 
% The cost equation can be estimated. The default intercept in generated data is 4 and slope
% (the |gamma| variable) has been set to 1 in the creation of the dataset.

market = Market(demand);
market.var.firm = 'firm';
market.settings.conduct = 0.5;

market.findCosts();
market.y = market.c;
market.var.exog = 'w';
market.var.panel = 'productid';
costEstimate = market.estimate()