%% SimMarket Users Manual

%% Demand estimation
% Demand is estimated using the |NestedLogitDemand| and
% |MixedLogitDemand| classes. The classes have methods to estimate demand of
% the respective type. 
%
% Data used in estimation is contained in a Matlab table object. Tables
% in Matlab are very similar to datasets in Stata or data frames in R. They
% can contain categorical/factor and numerical variables. Reference to
% variables in the table are by their variable names. 
%
% Demand is estimated using the |NestedLogitDemand| and |MixedLogitDemand|
% classes.

%% NestedLogitDemand estimation
% We will start by describing estimation of nested logit demand. Although the
% syntax is a little different than with _mergersim_ in Stata, the required fields are
% the same. We load a data table |dt|, and provide this table as an input to
% the demand constructor:

load example_data;
demand = NestedLogitDemand(dt1);

%%
% Here we specify the market and panel id, price, quantity and marketsize
% (here set to the constant value 1). In addition one can specify a list of
% exogenous variables to be used in estimation in |demand.var.exog|. 

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';

%%
% Properties can be defined in any order. To specify several exogenous
% variables, specify them in a string separated by spaces.
%
% There are a number of parameters that can be set in |NestedLogitDemand|
% The most important characteristics are |demand.var| and |demand.settings|. In
% |demand.var|, various variable names in the dataset are specified.

demand.var

%%
% In |demand.settings|, various settings are set. Four properties
% of |NestedLogitDemand|.settings concern estimation. The last one,
% |demand.settings.ces| is used to select CES Demand rather than the default Unit demand.

demand.settings

%%
% To estimate demand, the |estimate()| method is used. The estimation method 
% used will depend on what has been specified in demand.var and in
% demand.settings. 

result = demand.estimate()

%%
% The method returns a table with the estimate, as well as putting it and
% various other information in a demand.results structure.

%% MixedLogitDemand estimation
% Estimation of mixed logit is rather similar to nested logit. There are
% more parameters that one can specify, however

load example_data;
demand = MixedLogitDemand(dt2);

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
% The demand class |MixedLogitDemand| has more settings than
% |NestedLogitDemand|

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

market = Market(demand);
market.var.firm = 'productid';

%%
% Here we specify single product firms (corresponding to the product id).
% To find the marginal costs that correspond with the estimated demand and
% prices and quantities, we use the function |findCosts()|:

market.findCosts();
averageCosts = mean(market.c)

%%
% Having determined costs, one can use the market class (with its
% associated demand) to study variations in ownership, costs etc. The
% simplest way to do this is to make a copy of the |Market| object |market| and
% compute a new equilibrium with the copy, |market2|. The effects of the change in
% the market conditions in the two settings can then be compared using |market.compare()|.

market2 = copy(market);
market2.firm(market2.firm == 2 ) = 1;
market2.equilibrium();
mergerResult = market.compare(market2.p)

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

m = SimMarket()

%%
% The |SimMarket| object m contains a structure of settings m.model. In
% addition to the associated demand object it creates a new demand object *m.estDemand*
% that is used for estimation. The dataset created by |SimMarket| is stored
% in m.data.

%% 
% A demand model can be created by using one of the classes
% |NestedLogitDemand| or |MixedLogitDemand|. The demand object is created as
% follows:

demand = NestedLogitDemand();
demand.alpha = 1;

%%
% This command creates an unnested logit demand object, as the only
% property set is demand.alpha.
%
% To create a simulated dataset with 100 observations based on the demand
% object, a |SimMarket| object is created, and the demand object is attached

%%
% To create the dataset the method |init()| can be used. Invoking
% this method, changes the |SimMarket| object we have created. The object
% |m| now contains a dataset |m.data|.

m.demand = demand;
m.init()
m

%%
% By default |m.init()| ceates a market with 5 products and 100 markets, in long 
% format as a Matlab table. 
% The first 10 observations of the dataset are shown below. It contains a 
% market and product identifiers, a constant, costs |c|, a demand characteristic |x| and
% a variable |d| containing both observable and unobservable characteristics. By
% default, the disturbances containing both an individual and a product
% specific shock, the latter uncorrelated with observables (random
% effects).

m.data(1:10,:)

%%
% To add random prices and the corresponding demand the method |calculateDemand()| 
% is used. Using the demand specification in |m.demand|, quantities are calculated 
% based on prices and product characteristics |d|. The total
% share of the product including the outside good is shown below, as are
% the average prices and quantities by product. (By default market size is
% set to 1 in all markets).

sresults = m.calculateDemand()

%%
% Demand can be estimated using the estimate() method. With the standard
% configuration used here, this gives an OLS panel estimate, using fixed
% effects. The results shown have true values in the first column and the
% estimated values in the other columns. 

results = m.estimate()

%% MixedLogitDemand Monte-Carlo
% Now we will create a slightly more complex market with mixed logit demand. 
% A minimal definition of a |MixedLogitDemand| object is as follows:

demand = MixedLogitDemand();
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
% Here we let prices be endogenous. As an alternative to |calculateDemand()| 
% presented above, we instead use |simulateDemand()|. Instead of prices being
% random, |simulateDemand()| calculates equilibrium values depending on market conditions in each
% market. Price variability can comes from cost shifters and/or the number of products being set to be 
% exogenously random. Prices and quantities will depend on the products in
% the market as well as the ownership structure. 

m2.model.endog = true;
m2.model.randproducts = true;
m2.init();
m2.simulateDemand()

%%
% A minimal nested logit demand specification is as follows:

demand = NestedLogitDemand();
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

% By default SimMarket assumes single product firms. We can change this mapping
% by setting the |model.firms| property. Here we assume that the five
% products in the model have two owners. 

m3.model.firm = [1,1,1,2,2]';

m3.model.endog = true;
m3.model.randproducts = true;

m3.model.gamma = 1;
m3.init();
m3.simulateDemand()
dt3 = m3.data;

%%
% To estimate this model, the nesting variable |type| has to be specified.
% The same count instruments as above are used.

demand = NestedLogitDemand(dt3);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';

demand.var.nests = 'type';
demand.var.instruments = 'nprod nprod2';
demand.settings.paneltype = 'lsdv';
result = demand.estimate()

%%
% The datasets that have been created can be saved for later use, here to 
% the file |example_data.mat|:

dt1 = m.data;
dt2 = m2.data;
save example_data dt1 dt2 dt3;

