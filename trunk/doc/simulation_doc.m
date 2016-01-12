%% SimMarket 0.2 Monte-Carlo Simulation 

%% SimMarket class
% To create a simulated market, the |SimMarket| class is used. It has methods
% to create price and quantity variables:
% # from a random price using a specified demand
% # from a set of products that vary exogenously over time
% 
% In both cases instruments are created.
% 
% The purpose of the |SimMarket| class is to 
% # create a dataset based on a demand model and to 
% # facilitate estimation by creating a new demand object associated with the data. 

m1 = SimMarket()

%%
% The |SimMarket| object |m1| contains a structure of settings |m1.model|. 
% The dataset created by |SimMarket| is stored in |m.data|.

%% 
% A demand model can be created by using one of the classes
% |NLDemand| or |RCDemand|. The demand object is created as
% follows:

demand = NLDemand();
demand.alpha = 1;

%%
% This command creates an unnested logit demand object, as the only
% property set is demand.alpha.
%
% To create a simulated dataset with 100 observations based on the demand
% object, a |SimMarket| object is created, and the demand object is attached

m1.demand = demand;

%%
% To create the dataset the method |create()| is used. Invoking
% this method, changes the |SimMarket| object we have created. 

sresults = m1.create()

%%
% The object |m| now contains a dataset |m1.data|.
% By default |m1.create()| ceates a market with 5 products and 100 markets, in long 
% format as a Matlab table. 
% The first 10 observations of the dataset are shown below. It contains a 
% market and product identifiers, a constant, costs |c|, a demand characteristic |x| and
% a variable |d| containing both observable and unobservable characteristics. By
% default, the disturbances containing both an individual and a product
% specific shock, the latter uncorrelated with observables (random
% effects).

%%
% Using the demand specification in |m1.demand|, quantities are calculated 
% based on prices and product characteristics |d|. The total
% share of the product including the outside good is shown below, as are
% the average prices and quantities by product. (By default market size is
% set to 1 in all markets).

m1.data(1:10,:)

%%
% Once |m1.create()| has been run, a data table has been created for the
% market that can be used in estimation.

dt1 = m1.data;

%%
% |SimMarket| has an |estimate()| method that can be used to estimate
% demand based on the data created. It estimates the demand based on a copy
% of the demand object that has specified. 
%
% Demand can be estimated using the estimate() method. With the standard
% configuration used here, this gives an OLS panel estimate, using fixed
% effects. The results shown have true values in the first column and the
% estimated values in the other columns. This is useful in testing
% estimation methods.

results = m1.estimate()

%% RCDemand Monte-Carlo
% Now we will create a slightly more complex market with mixed logit demand. 
% A minimal definition of a |RCDemand| object is as follows:

demand = RCDemand();
demand.alpha = 1;
demand.sigma = 1;
demand.var.nonlinear = 'x';

%% 
% The |SimMarket| class has a set of model settings that can be used to
% customize the simulated market. All model settings have default values:

m2 = SimMarket();
m2.demand = demand;
display(m2.model)

%% 
% In contrast with the previous example, here we let prices be endogenous. 
% Instead of prices being exogenously random, |create()| calculates 
% equilibrium values depending on market conditions in each
% market. Price variability can comes from cost shifters and/or the number of products being set to be 
% exogenously random. Prices and quantities will depend on the products in
% the market as well as the ownership structure. Here we set ownership of
% the five products explicitly with the |m2.model.firm| setting. Note that
% as |m2.model.randomProducts = true| is specified, not all five products
% will actually exist in all markets.

m2.model.endog = true;
m2.model.randomProducts = true;
m2.model.firm = [1,1,2,2,3];
m2.create();

%%
% The data table created is as above stored in |m2.data|. 

dt2 = m2.data;

%% Nested logit Monte-Carlo
% A minimal nested logit demand specification with the nesting variable
% |type| is created as follows:

demand = NLDemand();
demand.alpha = 0.5;
demand.sigma = 0.5;
% demand.var.nests = 'type';

%%
% The market created can be specified by modifying the model characteristics. 
% To allow for nests, we can add a categorical variable |type| 
% having 2 distinct values by specifying:

m3 = SimMarket();
m3.demand = demand;
m3.model.types = 2;

%%
% Prices in |SimMarket| are by default determined by
% equilibrium conditions, with variation coming from changes in costs
% and/or the number of products in the market. 
%
% Alternatively, one can let prices be set randomly, with or without
% correlation with the error term. In this setup, quantities are calculated
% based on prices that are randomly generated. With endogeneity, a set of
% valid instruments inst1-inst6 are also generated.

%  m3.model.pricesFromCosts = false; 

%%
% By default SimMarket assumes single product firms. We can change this mapping
% by setting the |model.firms| property. Here we assume that the five
% products in the model have two owners. 

m3.model.firm = [1,1,1,2,2];

m3.model.endog = true;
% m3.model.randomProducts = true;

%%
% We add a cost shifter |w| by specifying the |model.gamma| parameter.

m3.model.gamma = 1;
m3.create();
dt3 = m3.data;

%%
% To estimate this model, the nesting variable |type| has to be specified.
% The same count instruments as above are used. Note that the 2SLS FE panel
% estimate will have price and log group shares as endogenous variables.

demand = NLDemand(dt3);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';

demand.var.nests = 'type';
demand.var.instruments = 'c w';
demand.settings.paneltype = 'lsdv';
result = demand.estimate()

%% 
% The cost equation can be estimated. The default intercept in generated data is 4 and slope
% (the |gamma| variable) has been set to 1 in the creation of the dataset.

market = Market(demand);
market.var.firm = 'firm';
market.findCosts();
market.y = market.c;
market.var.exog = 'w';
market.var.panel = 'productid';
costEstimate = market.estimate()

%%
% The datasets that have been created can be saved for later use, here to 
% the file |example_data.mat|:

save example_data dt1 dt2 dt3;

