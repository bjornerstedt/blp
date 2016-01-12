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
demand2 = RCDemand(dt2);

demand2.var.market = 'marketid';
demand2.var.panel = 'productid';
demand2.var.price = 'p';
demand2.var.quantity = 'q';
demand2.var.marketsize = 'constant';
demand2.var.exog = 'x';
demand2.var.instruments = 'nprod nprod2';

demand2.var.nonlinear = 'x';

%%
% To estimate, the |estimate()| method is invoked:

result = demand2.estimate()

%%
% The demand2 class |RCDemand| has more settings than
% |NLDemand|

demand2.settings

%% 
% We can for example estimate using quadrature, with optimal instruments:

demand2.settings.drawmethod = 'quadrature';
demand2.settings.optimalIV = true;
result = demand2.estimate()

%% Market class
% The market class is used to calculate equilibrium prices and quantities,
% based on market structure. The equilibrium depends on the estimated
% parameters of the demand model specified.

market = Market(demand2);
market.var.firm = 'productid';

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
