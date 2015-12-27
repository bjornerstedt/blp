%% SimMarket Reference

%% Estimate class
% SimMarket demand and market classes all inherit the linear estimation 
% functionality of the |Estimate| class. This class can be used for
% estimation not directly related to demand or market estimation. 

load example_data;
est = Estimate(dt3)
est.var

%%
%   Estimate with properties:
% 
%         data: []
%      panelid: []
%     marketid: []
%      results: [1x1 struct]
%            y: []
%            X: []
%            Z: []
%         beta: []
%     settings: [1x1 SettingsClass]
%       config: [1x1 struct]
%          var: [1x1 SettingsClass]
%        Xorig: [x]
%        Zorig: [x]
%        
% est.var:
%          market: []
%           panel: []
%          depvar: []
%            exog: []
%           endog: []
%     instruments: []
    
est.settings

%%
% est.settings
%             robust: 'true'
%          paneltype: 'fe'/'lsdv'/'none'
%             nocons: 0
%            weights: []
%     estimateMethod: 'ols'/'2sls'/'gmm'
%
% est.results
%     estimateDescription: 'Linear Estimate'
%                   other: [x]
%                  params: [1x1 struct]
%                estimate: [2x3 table]
%                     var: [6x2 table]
%                settings: [5x2 table]
%                
%     est.settings.paneltype = 'none';
%     est.var.exog = 'w';
%     est.var.depvar = 'c';
%     est.estimate()

%% NestedLogitDemand class
%

demand = NestedLogitDemand(dt1)
demand.var

%%
%   NestedLogitDemand with properties:
% 
%        alpha: []
%        sigma: []
%            d: []
%         data: [500x9 table]
% 
%   demand.var:
% 
%           price: []
%           nests: []
%        quantity: []
%      marketsize: []
%% 
% Settings
% 
%                ces: 0

demand.settings

%% 
% Methods

%% MixedLogitDemand class
%

demand = MixedLogitDemand(dt1)
demand.var

%%
%  MixedLogitDemand with properties:
% 
%     rc_sigma: []
%        draws: []
%           x2: [x]
%            W: []
%           xi: [x]
%        
%% 
% Settings
% 
%   SettingsClass with properties:
% 
%            maxiter: 100
%       fptolerance1: 1.0000e-14
%          quaddraws: 10
%          optimalIV: 0
%        marketdraws: 0
%       fptolerance2: 1.0000e-14
%         drawmethod: 'hypercube'
%                ces: 0
%               nind: 100

demand.settings

%%
% demand.config
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
% Methods


%% Market class
% 
%         firm: []
%            q: []
%            p: []
%           p0: []
%            c: []
           
market = Market()
market.var

%% 
% Settings
% 
%               dampen: 1
%              conduct: 0
%     weightedAverages: 1
%              weights: []
%          valueShares: 0

market.settings

%% 
% Methods


%% SimMarket class
% In addition to the associated demand object it creates a new demand object |m.estDemand|
% that is used for estimation. 

m = SimMarket()


%% 
% m.model
%
%              endog: 0
%       randproducts: 0
%     simulatePrices: 1
%            markets: 100
%           products: 5
%              types: []
%               firm: []
%               beta: [1 0]
%                  x: [5 0]
%            x_sigma: [1 1]
%                  c: 4
%            c_sigma: 1
%              gamma: 0
%      epsilon_sigma: 0.1000
%           sigma_xi: 0.1000
%        endog_sigma: 0.1000
%          prob_prod: 0.8000
         
m.model

%% 
% Methods

