%function bethorg = createmarket(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP ESTIMATION WITH SIMULATED DATA SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ;
% if nargin == 0
     sigmaEpsilon = 0.1
% else
%     sigmaEpsilon = varargin{1};
% end
% sigmaEpsilon

randproducts = true 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATION OF DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('*****  Generating demand  *****');
%SET SEED
stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

% INITIAL SETTING (SIMILAR TO COFFEE DATA)
nmkt = 500;   % number of markets = (# of cities)*(# of quarters)  
nbrn = 8;     % number of brands per market

nevo = table();
n = nmkt*nbrn;  % number of observations    
nevo.cdid = kron((1:nmkt)',ones(nbrn,1)); %vector with marketnumber for each obs
nevo.brn = mod(0:(n - 1), nbrn)' + 1; 

% DEFINE TRUE PARAMETERS
betatrue = [ -1; 1; -1.8; 30];
%sigmaEpsilon = 0.05;           %standard deviation of individual prefs
rc_sigmatrue = [1];


xi = randn(n,1) * sigmaEpsilon;  %unobserved characteristics
nevo.f = 12 + randn(n, 1);     % feature
nevo.d = 5 + randn(n, 1);     % display

%Simulate Instruments
A = 2 + 0.2 * randn(n, 6);
M = ones(6, 6)*0.8;
M(1:7:36) = 1; % Make diagonal elements = 1
instruments = A*chol(M);

%Simulate endogenous price. Note f and d do not affect price.
nevo.p = 5.5 + xi + sum(instruments, 2); % Endogenous price
%nevo.p = 5.5 + randn(n,1) + sum(instruments, 2); % Exogenous price
%nevo.p = 5.5 + sum(instruments, 2); % Exogenous price
%Brand Dummies
%brand=repmat([1:nbrn]', nmkt, 1);
%branddum=dummyvar(brand);
nevo.ms = ones(n, 1);
if randproducts
    % prob .8 of selecting product in a period
    prodsel = logical(binornd(1,.8, n, 1));
    nevo = nevo(prodsel, :);
    xi = xi(prodsel, :);
end
% Create count instrument
nprod = accumarray(nevo.cdid, nevo.ms);
nprod = nprod(nevo.cdid,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION OF DEMAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

demand = MixedLogitDemand(nevo);
demand.var,exog = 'f d';
% demand.endog = [ ];
demand.var.panel = 'brn';
demand.var.market = 'cdid';
demand.p = 'p';
demand.var.price = 'p';
demand.beta = betatrue;
demand.rc_sigma = rc_sigmatrue;
demand.xi = xi;
demand.ispanel = false;
demand.halton = true;
demand.halton2 = false;
demand.quadrature = false;
demand.nind = 100;
demand.marketdraws = false;
demand.fptolerance1 = 1e-12;
demand.quaddraws = 15;
demand.nonlinear = 'f';

demand.init();
draws = demand.draws;
save draws draws;

display('*****  Generating market  *****');
market = Market(demand);
%market.firm = nevo.brn;
market.firmvar = 'brn'; % One product firms

if false
    market.s = demand.s; 
    market.p0 = demand.p;
    market.p = demand.p;
    market.init();
    market.findCosts();
    res = [market.c, market.p]
    mean(res,1)
    var(res,1)
else
    cmean = 20;
    cvar = 2;
    market.c = cmean + cvar*randn(size(market.marketid));
%     market.p0 = demand.p;
%     market.p = demand.p;
    market.p0 = market.c;
    market.p = market.c;
    market.init();
    found = market.equilibrium();
end
res = [market.c, market.p];
res
nevo.shares  = demand.shares(nevo.p);
nevo.p = market.p;
bethorg = [betatrue ; rc_sigmatrue];

inst = [market.c, nprod, nprod .^2];
nevo = [nevo array2table(inst)];
nevo = nevo(found(nevo.cdid,2) > 0, : );

save simmarket nevo bethorg xi;
display('*****  Finished  *****');

