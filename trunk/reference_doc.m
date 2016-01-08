%% SimMarket 0.2 Reference 
% This is a preliminary version that includes output of structures for
% comparisons.

%% Estimate class
% SimMarket demand and market classes all inherit the linear estimation 
% functionality of the |Estimate| class. This class can be used for
% estimation not directly related to demand or market estimation. 
load example_data;
est = Estimate(dt3);

%%
% The Estimate class has the following properties
% 
%    settings: A structure with different estimation settings
%     config: Contains less common settings
%        var: A structure with variable names used in estimation
%       data: The Matlab table with data useed in estimation. Can be
%             specified in the constructor as above.
%    results: A structure with results (coefficients, standare errors, other statistics) 
est

%%        
% est.var contains fields for variables used in estimation:
% 
%          market: Misnomer for this general class - change to time?
%           panel: Panel data identifier
%          depvar: Dependent variable
%            exog: List of exogenous variable names, separated by spaces
%           endog: List of endogenous variable
%     instruments: List of instruments   
est.var

%%
% The est.settings structure has the following fields
%
%             robust: 1 - robust estimation true/false
%          paneltype: 'none' - panel estimate: 'fe'/'lsdv'/'none'
%             nocons: 0 Do not include constant in estimation true/false
%     estimateMethod: 'ols'/'2sls'/'gmm'
est.settings

%% 
% The method |Estimate.estimate()| generates a result table as output. It
% also populates the structure |Estimate.results| with various results
%
% estimateDescription: 'Linear Estimate'
%               other: [x]
%              params: Structure with the estimate and var-covar matrix
%            estimate: Estimate table
%                 var: Table with variable names used
%            settings: Table with settings
est.settings.paneltype = 'none';
est.var.exog = 'w';
est.var.depvar = 'c';
est.estimate()

%%
% The complete variance covariance matrix is obtained in the property
est.results.params.varcovar
%%
% The settings and var structs can also be obtained as tables
est.results.var
est.results.settings

%%
% In estimating, Estimate creates the following properties that can be used in
% calculations. To generate these properties without estimating, the method
% |Estimate.init()| is invoked.
% 
%         y: []
%         X: []
%         Z: []
%      beta: []
%
%   panelid: []
%  marketid: []
%     Xorig: [x]
%     Zorig: [x]


%%
% *Methods*
%
% Estimation is done with the |estimate()| method. The mehod used depends
% on the type of object that estimation is performed on. In the |Estimate|
% class, the method can be set to OLS, 2SLS or GMM in settings. 
methods(Estimate)

%% NestedLogitDemand class
%
% The demand classes extend |Estimate| to allow estimation of demand
% systems. 
demand = NestedLogitDemand(dt1)
demand.var

%%
% NestedLogitDemand has the following additional properties:
% 
%      alpha: The calibrated or estimated alpha parameter
%      sigma: A vector with sigmas
%          d: A vector with utility shifters, used in Monte Carlo estimation
% 
% Additional variables are specified in demand.var:
% 
%         price: Variable name of price variable
%         nests: Name(s) of nesting variables
%      quantity: Quantity variable
%    marketsize: Name of variable in dataset containing market size per market
%
%
% There is also an additional setting in NestedLogitDemand beyond those of
% Estimate:
% 
% * ces: 0 - Use CES logit rather than unit demand true/false
demand.settings

%% 
% *Methods*
%
% The method |NestedLogitDemand.estimate()| performs a linear panel
% estimate based on the settings.
methods(NestedLogitDemand)

%% MixedLogitDemand class
%
demand = MixedLogitDemand(dt1)
demand.var

%%
% MixedLogitDemand with properties:
% 
%    rc_sigma: The calibrated or estimated nonlinear parameters
%        
%% 
% Settings
% 
% MixedLogitDemand.settings has properties:
% 
%              ces: 0 - CES or Unit logit demand
%          maxiter: 100 - Maximum number of iterations in optimization
%        optimalIV: 0 - Optimal instruments true/false
%       drawmethod: 'hypercube' - Sampling method:
%       'hypercube'/'quadrature'/'halton'/'random'
%             nind: 100 - Number of simulated individuals
%      marketdraws: 0 - Different random draws for each market true/false
%        quaddraws: 10 - Quadrature accuracy level
%     fptolerance1: 1.0000e-14
%     fptolerance2: 1.0000e-14
demand.settings

%%
% MixedLogitDemand.config
% 
%                  hessian: 0
%                     test: []
%                  fpmaxit: 1000
%                tolerance: 1.0000e-09
%               randstream: []
%              restartFval: 1000
%               guessdelta: 1
%                  quietly: 1
%     restartMaxIterations: 1
demand.config

%% 
% *Methods*
%
% The method |MixedLogitDemand.estimate()| performs a BLP
% estimate based on the settings specified in the demand object.
methods(MixedLogitDemand)

%% Market class
% 
% The |Market| class is used to calculate costs or to 
% associated with a demand class either in its
% constructor or by setting |Market.demand|
%
% demand: Demand object (|NestedLogitDemand| or |MixedLogitDemand|)
%      p: Equilibrium calculated price
%      q: Equilibrium calculated quantity
%     p0: Initial guess for equilibrium price
%      c: Costs calculated from market prices and quantities and demand estimate
% 
% The |Market| class obtains data and various settings from the associated
% demand class. List these... It has the settings and var structures allowing estimation of costs. 
market = Market();
market.var

%%
% |Market.var.firm| is the only property that has to be set to calculate
% costs or equilibrium.
market.var.firm = 'productid';
%%
% An estimated or calibrated demand is associated with the marktet either
% by providing it in the constructor or adding it to the |Market.demand|
% property:
market.demand = demand;

%% 
% The Market class has the following settings, set in Market.settings
% 
%           dampen: 1 - Dampening in fixed point iterations
%            maxit: 1000 - Maximum number of iterations in calculating equilibrium
%          conduct: 0 - Conduct parameter in [0,1] interval
% weightedAverages: 1 - Calculate weighted averages true/false
%      valueShares: 0 (1 for CES) - Use value shares as weights true/false
market.settings

%% 
% *Methods*
%
% |Market.findCosts()| calculates costs based on a demand specification
% Prices and quantities used are copied from the demand specification
% 
%   market.findCosts()

%%
% |Market.equilibrium()| calculates a market equilibrium based on a demand
% specification, costs, and a specification of ownership and conduct (using
% |Market.var.firm| and |Market.settings.conduct|. 
methods(Market)

%%
% To summarise market calculations after calculating costs with |Market.findCosts()|, 
% |Market.summarise()| can be used. By default weighted averages by firm using market
% shares as weights are calculated. Market shares by quantity are used by
% default for unit demand, and market shares by value for CES demand. The
% type of average can be controlled by setting
% |Market.settings.weightedAverages| and |Market.settings.valueShares|.
% One can also summarise a selection in a logical vector by invoking:
%
%   m1.summary('selection', vec);
% 
% Options to functions in Matlab are specified with argument pairs: an
% identifying text string and a value. To get averages by product for the
% same selection, one specifies:
%
%   m1.summary('selection', vec, 'GroupingVariables', 'productid'); 

%%
% The method |Market.compare()| can be used to compare two market equilibria 
% m1 and m2. Aggregated prices and percentage price changes from m1 to m2 are shown. 
% By default weighted average prices, grouped by firm are shown using
%
%   compare(m1, m2)
% 
% |compare| can take the same options as |summarise|.

%% SimMarket class
% In addition to the associated demand object it creates a new demand object |m.estDemand|
% that is used for estimation. 
%
% SimMarket has the following properties:
% 
%  model: Structure with model settings
%   data: Data created by SimMarket
% demand: Demand model specified by user
% market: Market model specified by user
m = SimMarket()

%% 
% The market model is specified with the parameters in |SimMarket.model|:
% 
%          endog: 0 - Endogenous prices and quanities true/false
%   randproducts: 0 - Exogenously random products in market true/false
% simulatePrices: 1 - Simulate prices or let them be randomly drawn as
%                     in Nevo code true/false
%        markets: 100 - Number of markets generated
%       products: 5 - (Maximum) number of products in each market.
%          types: [] - Number of types for each categorical
%           firm: [] - Vector of ownership for each producty
%           beta: [1 0] - 
%              x: [5 0] - Expected value for p and other demand shifters
%          x_vcv: [1 1] - Variance, can be specified as a matrix for
%                         multicollinearity
%              c: 4 - Costs
%          c_vcv: 1
%          gamma: 0 - Cost shifter parameter
%  epsilon_sigma: 0.1 - Sd of individual unobservables
%       sigma_xi: 0.1 - Sd of panel unobservables
%    endog_sigma: 0.1 - Endogeneity parameter for non simulated prices
%      prob_prod: 0.8 - Probability that product exists in a market       
m.model

%% 
% *Methods*
% 
% * SimMarket - Create a new simulation object, optionally with demand spec       
% * create - Creates market - should return dataset.
% * estimate -  Estimate and compare, used in testing framework
% * findCosts - Calculate costs, used in testing framework
methods(SimMarket)

