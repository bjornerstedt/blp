%function bethorg = create(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP ESTIMATION WITH SIMULATED DATA SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETERMINATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear ;
%     if nargin == 0
        sigmaEpsilon = 0.8
        sigmaEpsilon = 0.05
%     else
%         sigmaEpsilon = varargin{1};
%     end

mixed = false

display('*****  Generating demand  *****');
for ces = 0:1

stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

% INITIAL SETTING (SIMILAR TO COFFEE DATA)
nmkt = 400;   % number of markets = (# of cities)*(# of quarters)  
nbrn = 6;     % number of brands per market

nevo = table();
n = nmkt*nbrn;  % number of observations    
nevo.cdid = kron((1:nmkt)',ones(nbrn,1)); %vector with marketnumber for each obs
nevo.brn = mod(0:(n - 1), nbrn)' + 1; 

% DEFINE TRUE PARAMETERS
betatrue = [ -1; 1; -1.8; 1.9];
%sigmaEpsilon = 0.05;           %standard deviation of individual prefs
rc_sigmatrue = [1];
%rc_sigmatrue = [1; 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATION OF DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi = randn(n,1) * sigmaEpsilon;  %unobserved characteristics
nevo.f = 12 + randn(n, 1);     % feature
nevo.d = 5 + randn(n, 1);     % display

%Simulate Instruments
A = 2 + 0.2 * randn(n, 6);
M = ones(6, 6)*0.8;
M(1:7:36) = 1; % Make diagonal elements = 1
inst = A*chol(M);

%Simulate endogenous price. Note f and d do not affect price.
nevo.p = 5.5 + xi + sum(inst, 2); % Endogenous price
%nevo.p = 5.5 + randn(n,1) + sum(instruments, 2); % Exogenous price
%nevo.p = 5.5 + sum(instruments, 2); % Exogenous price
%Brand Dummies
%brand=repmat([1:nbrn]', nmkt, 1);
%branddum=dummyvar(brand);
constant=ones(n, 1);
nevo = [nevo array2table(inst)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION OF DEMAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mixed
    if ces == 1
        demand = CesMixedLogitDemand(nevo);
    else
        demand = MixedLogitDemand(nevo);
    end
    demand.beta = betatrue;
    demand.rc_sigma = rc_sigmatrue;
else
    if ces == 1
        demand = CesNestedLogitDemand(nevo);
    else
        demand = NestedLogitDemand(nevo);
    end
    demand.alpha = -2;
end
demand.exog = 'f d';
demand.endog = [ ];
demand.panel = 'brn';
demand.market = 'cdid';
demand.p = 'p';
demand.price = 'p';
demand.xi = xi;
demand.settings.paneltype = 'none';
demand.halton = true;
demand.quadrature = false;
demand.nind = 300;
demand.marketdraws = false;

demand.nonlinear = 'f';

demand.init();
draws = demand.draws;
save draws draws;
bethorg = [betatrue ; rc_sigmatrue];

nevo.shares  = demand.shares(nevo.p);
if ces == 1
    MS= 100;
nevo.ms = constant .* MS;
    nevo.quantity = nevo.shares ./ nevo.p * MS;
    save simdatacesX nevo bethorg xi ces;
else
nevo.ms = constant;
    nevo.quantity  = nevo.shares;
    save simdataX nevo bethorg xi ces;
end
%{
deltatrue = demand.x1*betatrue + xi;
if size(demand.x2,1) > 0,
    chamb = numgrad(@demand.delta, rc_sigmatrue); % ??
end
instrument = [instrument chamb];
%}
end
display('*****  Finished  *****');

