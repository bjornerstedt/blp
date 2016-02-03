
%% Test 1: Create data
% Test creation of simulated data
% Invoke with: runtests simtest1

ces = false

display('*****  Generating demand  *****');

stream1 = RandStream('mt19937ar', 'Seed', 999);
RandStream.setGlobalStream(stream1);

nmkt = 250;   % number of markets = (# of cities)*(# of quarters)  
nbrn = 16;     % number of brands per market

nevo = table();
n = nmkt*nbrn;  % number of observations    
nevo.marketid = kron((1:nmkt)',ones(nbrn,1)); %vector with marketnumber for each obs

% DEFINE TRUE PARAMETERS
betatrue = [ -1; 1; 1];
sigmaEpsilon = 0.001;           %standard deviation of individual prefs
rc_sigmatrue = [1];

xi = randn(n,1) * sigmaEpsilon;  %unobserved characteristics
% nevo.f = 12 + randn(n, 1);     % feature
nevo.d = 5 + randn(n, 1);     % display
nevo.p = 5 + randn(n, 1);
%Brand Dummies
%brand=repmat([1:nbrn]', nmkt, 1);
%branddum=dummyvar(brand);
constant=ones(n, 1);
nevo.ms = constant;

demand = MixedLogitDemand(nevo);
demand.settings.ces = ces;
demand.settings.marketdraws = true;
demand.var.marketsize = 'ms';
demand.var.price = 'p';
demand.var.nonlinear = 'd'; 

demand.var.market = 'marketid';
%demand.var.panel = 'product';
demand.var.exog = 'd';
demand.settings.paneltype = 'none';

demand.settings.nind = 1000;
demand.settings.drawmethod = 'hypercube';
% demand.settings.drawmethod = 'quadrature';
% obj.settings.quaddraws = 12;

demand.beta = betatrue;
demand.rc_sigma = rc_sigmatrue;
demand.init(); 

demand.initSimulation();

% draws = demand.draws;
% save draws draws;
% bethorg = [betatrue ; rc_sigmatrue];

nevo.shares  = demand.shares(nevo.p, xi);
assert(~any(isnan(nevo.shares)))

% Test 2: Estimate on simulated data
% Test both that results do not change, and that close to true value

display('*****  Estimating  *****');
demand = MixedLogitDemand(nevo);
demand.settings.marketdraws = true;

demand.var.quantity = 'shares';
demand.settings.ces = ces;
demand.var.marketsize = 'ms';
demand.var.price = 'p';
demand.var.instruments = '';
demand.var.nonlinear = 'd'; 

demand.var.market = 'marketid';
%demand.var.panel = 'product';
demand.var.exog = 'd';
demand.settings.paneltype = 'none';

demand.settings.nind = 1000;
demand.settings.drawmethod = 'hypercube';
% demand.settings.drawmethod = 'quadrature';
% obj.settings.quaddraws = 12;
demand.config.tolerance = 10e-12;
obj.settings.fptolerance1 = 1e-16; % lower tol for first FP its
obj.settings.fptolerance2 = 1e-16; % high tol for last iterations

            
demand.init(); 
demand.estimate(); 

assert(~isempty(demand.results.estimate))

