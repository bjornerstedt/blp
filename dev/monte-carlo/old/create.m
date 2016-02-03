%function bethorg = create(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE SIMULATED DATA SET FOR BLP ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ;
%     if nargin == 0
%         sigmaEpsilon = 0.8
%     else
%         sigmaEpsilon = varargin{1};
%     end

ces = true

display('*****  Generating demand  *****');

stream1 = RandStream('mt19937ar', 'Seed', 999);
RandStream.setGlobalStream(stream1);

% INITIAL SETTING (SIMILAR TO COFFEE DATA)
nmkt = 200;   % number of markets = (# of cities)*(# of quarters)  
nbrn = 4;     % number of brands per market

nevo = table();
n = nmkt*nbrn;  % number of observations    
nevo.marketid = kron((1:nmkt)',ones(nbrn,1)); %vector with marketnumber for each obs

% DEFINE TRUE PARAMETERS
betatrue = [ -.1; .1; .1];
betatrue = [ -1; 1; 1];
sigmaEpsilon = 0.1;           %standard deviation of individual prefs
rc_sigmatrue = [1];

%% GENERATION OF DATASET

xi = randn(n,1) * sigmaEpsilon;  %unobserved characteristics
% nevo.f = 12 + randn(n, 1);     % feature
nevo.d = 5 + randn(n, 1);     % display
%Simulate Instruments
A = 2 + 0.2 * randn(n, 6);
M = ones(6, 6)*0.8;
M(1:7:36) = 1; % Make diagonal elements = 1
inst = A*chol(M);

%Simulate endogenous price. Note f and d do not affect price.
nevo.p = 5.5 + xi + sum(inst, 2); % Endogenous price
nevo = [nevo array2table(inst)];

% nevo.p = 5 + randn(n, 1);
%Brand Dummies
%brand=repmat([1:nbrn]', nmkt, 1);
%branddum=dummyvar(brand);
constant=ones(n, 1);
nevo.ms = constant;

%% SIMULATION OF DEMAND

demand = MixedLogitDemand(nevo);
demand.settings.ces = ces;
demand.var.marketsize = 'ms';
demand.var.price = 'p';
demand.var.nonlinear = 'd'; 

demand.var.quantity = 'ms';
demand.var.market = 'marketid';
% demand.var.panel = 'product';
demand.var.exog = 'd';
demand.settings.paneltype = 'none';

demand.settings.nind = 300;
demand.settings.drawmethod = 'hypercube';

demand.beta = betatrue;
demand.rc_sigma = rc_sigmatrue;
demand.init(); 

demand.initSimulation();

% draws = demand.draws;
% save draws draws;

% bethorg = [betatrue ; rc_sigmatrue];

nevo.shares  = demand.shares(nevo.p);
if ces == 1
    nevo.quantity = nevo.shares ./ nevo.p;
    save simdataces nevo betatrue xi ces;
else
    nevo.quantity  = nevo.shares;
    save simdata nevo betatrue xi ces;
end
%{
deltatrue = demand.x1*betatrue + xi;
if size(demand.x2,1) > 0,
    chamb = numgrad(@demand.delta, rc_sigmatrue); % ??
end
instrument = [instrument chamb];
%}
display('*****  Finished  *****');

