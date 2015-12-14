%% SimMarket Documentation

%% Demand estimation
% Demand is estimated using the NestedLogitDemand and
% MixedLogitDemand classes. The classes have methods to estimate demand of
% the respective type. 
%
% Data used in estimation is contained in a Matlab table object. Tables
% in Matlab are very similar to datasets in Stata or data frames in R. They
% can contain categorical/factor and numerical variables. Reference to
% variables in the table are by their variable names. 
%
% Demand is estimated using the NestedLogitDemand and MixedLogitDemand
% classes.

%% NestedLogitDemand estimation
% We will start by describing estimation of nested logit demand. Although the
% syntax is a little different than with *mergersim* in Stata, the required fields are
% the same. We load a data table dt, and provide this table as an input to
% the demand constructor:

load example_data;
demand = NestedLogitDemand(dt);

%%
% Here we specify the market and panel id, price, quantity and marketsize
% (here set to the constant value 1). In addition one can specify a list of
% exogenous variables to be used in estimation in demand.var.exog. 

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
% There are a number of parameters that can be set in NestedLogitDemand
% The most important characteristics are demand.var and demand.settings. In
% demand.var, various variable names in the dataset are specified.

demand.var

%%
% In demand.settings, various settings are set. The first four properties
% of NestedLogitDemand.settings concern estimation. The last one,
% settings.ces is used to select CES Demand rather than the default Unit demand.

demand.settings

%%
% To estimate demand, the estimate() method is used. The estimation method 
% used will depend on what has been specified in demand.var and in
% demand.settings. 

demand.estimate()

%%
% The method returns a table with the estimate, as well as putting it and
% various other information in a demand.results structure.

%% MixedLogitDemand estimation
% Estimation of mixed logit is rather similar to nested logit. There are
% more parameters that one can specify, however

load example_data;
demand = MixedLogitDemand(dt);

demand.var.market = 'marketid';
demand.var.panel = 'productid';
demand.var.price = 'p';
demand.var.quantity = 'q';
demand.var.marketsize = 'constant';
demand.var.exog = 'x';

demand.var.nonlinear = 'x';

%%
% To estimate, the estimate() method is invoked:

demand.estimate()

%%
% The demand class MixedLogitDemand has more settings than
% NestedLogitDemand

demand.settings

%% 
% We can for example estimate using quadrature, with optimal instruments:

demand.settings.drawmethod = 'quadrature';
% demand.settings.optimalIV = true;
demand.estimate()

%% Market class
% The market class is used to calculate equilibrium prices and quantities,
% based on market structure. The equilibrium depends on the estimated
% parameters of the demand model specified. 

market = Market(demand);
market.var.firm = 'productid';

%%
% Here we specify single product firms (corresponding to the product id).
% To find the marginal costs that correspond with the estimated demand and
% prices and quantities, we use the function findCosts():

market.findCosts();
display('Average costs are given by: ')
mean(market.c)

%% SimMarket Introduction
% A demand model can be created by using one of the classes
% NestedLogitDemand or MixedLogitDemand. The demand object is created as
% follows:

demand = NestedLogitDemand();
demand.alpha = 1;

%%
% This command creates an unnested logit demand object, as the only
% property set is demand.alpha.
%
% To create a simulated dataset with 100 observations based on the demand
% object, a SimMarket object is created, and the demand object is attached

m = SimMarket();

%%
% By default init() ceates a market with 5 products and 100 markets, in long 
% format as a Matlab table. 

%%
% To create the dataset the method init() can be used. Invoking
% this method, changes the SimMarket object we have created. The object m
% now contains a dataset.

m.demand = demand;
m.init()
m

%%
% The first 10 observations of the dataset are shown below. It contains a 
% market and product identifiers, a constant, costs c, a demand characteristic x and
% a variable d containing both observable and unobservable characteristics. By
% default, the disturbances containing both an individual and a product
% specific shock, the latter uncorrelated with observables (random
% effects).

m.data(1:10,:)

%%
% To add random prices and the corresponding demand the method calculateDemand() 
% is used. Using the demand specification in m.demand, quantities are calculated 
% based on prices and product characteristics d. The total
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


%% The MixedLogitDemand class
% A minimal definition of a MixedLogitDemand object is as follows:

demand = MixedLogitDemand();
demand.alpha = 1;
demand.rc_sigma = 1;
demand.var.nonlinear = 'constant';

%% The SimMarket class
% The SimMarket class has a set of model settings that can be used to
% customize the simulated market. All model settings have default values:

m = SimMarket();
m.demand = demand;
display(m.model)

%% 
% Here we let prices be endogenous. The number of products is set to be 
% exogenously random. Prices and quantities will depend on the products in
% the market as well as the ownership structure. By default we have single
% product firms

m.model.randproducts = true;
m.init();
m.simulateDemand()

%%
% The dataset used above is simply the data generated by the commands
% above:

dt = m.data;
save example_data dt;