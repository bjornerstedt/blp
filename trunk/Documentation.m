%% SimMarket 0.2 Users Guide

%% Demand estimation
% Demand is estimated using the |NLDemand| and
% |RCDemand| classes. The classes have methods to estimate demand of
% the respective type. 
%
% Data used in estimation is contained in a Matlab  <matlab:doc('table') table>  object. Tables
% in Matlab are very similar to datasets in Stata or data frames in R. They
% can contain categorical/factor and numerical variables. Reference to
% variables in the table are by their variable names. 
%
% Demand is estimated using the |NLDemand| and |RCDemand|
% classes.

%% NLDemand estimation
% We will start by describing estimation of nested logit demand. Although the
% syntax is a little different than with _mergersim_ in Stata, the required fields are
% the same. We load data including table |dt1|, and provide this table as an input to
% the demand constructor:

load example_data;
demand = NLDemand(dt1);

%%
% Here we specify the market and panel id, price, quantity and marketsize
% (here set to the constant value 1). In addition one can specify a list of
% exogenous variables to be used in estimation in |demand.var.exog|. A list
% of variables is specified as a string with variable names separated by
% spaces.

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';

%%
% There are a number of parameters that can be set in |NLDemand|
% The most important characteristics are set in |demand.var| and |demand.settings|. In
% |demand.var|, various variable names in the dataset are specified. One
% can list the contents of the structure:

demand.var

%%
% In |demand.settings|, other demand settings are set. Four properties
% of |NLDemand|.settings concern estimation. The last one,
% |demand.settings.ces| is used to select CES Demand rather than the default,
% Unit demand.

demand.settings

%%
% To estimate demand, the |estimate()| method is used. The estimation method 
% used will depend on what has been specified in demand.var and in
% demand.settings. 

result = demand.estimate()

%%
% The method returns a table with the estimate, as well as putting it and
% various other information in a demand.results structure.

%% RCDemand estimation
% Estimation of mixed logit is rather similar to nested logit. There are
% more parameters that one can specify, however. In the following example,
% we use count instruments to identify an endogenous price variable.

load example_data;
demand = RCDemand(dt2);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';
demand.var.instruments = 'nprod nprod2';

demand.var.nonlinear = 'x';

%%
% To estimate, the |estimate()| method is invoked:

result = demand.estimate()

%%
% The demand class |RCDemand| has more settings than
% |NLDemand|

demand.settings

%% 
% We can for example estimate using quadrature, with optimal instruments:

demand.settings.drawmethod = 'quadrature';
demand.settings.optimalIV = true;
result = demand.estimate()

%% Market class
% The market class is used to calculate equilibrium prices and quantities,
% based on market structure. The equilibrium depends on the estimated
% parameters of the demand model specified.

demand2 = copy(demand);
market = Market(demand);
market.var.firm = 'firm';

%%
% Here we specify single product firms (corresponding to the product id).
% To find the marginal costs that correspond with the estimated demand and
% prices and quantities, we use the function |findCosts()|:

market.findCosts();
averageCosts = mean(market.c)
market.summary()
market.summary('selection', dt2.marketid == 1)

%%
% Having determined costs, one can use the market class (with its
% associated demand) to study variations in ownership, costs etc. The
% simplest way to do this is to make a copy of the |Market| object |market| and
% compute a new equilibrium with the copy, |market2|. The effects of the change in
% the market conditions in the two settings can then be compared using |market.compare()|.

market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();
mergerResult = compare(market, market2)
mergerResult = compare(market, market2, 'selection', dt2.marketid == 1)

%%
% Cost calculation and equilibrium simulation can be performed on a
% selection rather than the whole dataset. To do this, a selection vector
% is provided. In this example we restrict our attention to market 1 by
% specifying |findCosts(dt2.marketid == 1)|. As costs have not been
% calculated for other markets, both |market2.equilibrium()| and |compare|
% calculate only for this market. 

market = Market(demand2);
market.var.firm = 'firm';
market.findCosts(dt2.marketid == 1);
market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();
mergerResult2 = compare(market, market2)

%%
% One can also explicitly restrict 
% equilibrium calculation to some markets. In this case calculations will
% only be on these markets, provided that costs have been calculated for
% them. Comparisons can be restricted to a subset of markets by providing
% the restriction with the |'Selection'| option.
%
%   market2.equilibrium(dt2.marketid == 1);
%   mergerResult2 = compare(market, market2, 'Selection', dt2.marketid == 1)

%%
% The |Market| class has a number of settings. Similarly to the demand
% classes, estimation is possible (see section below). One can also specify
% |market.settings.conduct| - the degree to which profits of other firms are taken in
% to consideration in profit maximization.

%% SimMarket Monte-Carlo
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
demand.rc_sigma = 1;
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
% as |m2.model.randproducts = true| is specified, not all five products
% will actually exist in all markets.

m2.model.endog = true;
m2.model.randproducts = true;
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
demand.var.nests = 'type';

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

m3.model.simulatePrices = false; 

%%
% By default SimMarket assumes single product firms. We can change this mapping
% by setting the |model.firms| property. Here we assume that the five
% products in the model have two owners. 

m3.model.firm = [1,1,1,2,2];

m3.model.endog = true;

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
demand.var.instruments = 'inst1 inst2 inst3 inst4 inst5 inst6';
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

