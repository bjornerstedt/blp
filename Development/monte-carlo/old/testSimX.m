%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLP ESTIMATION WITH SIMULATED DATA SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
    if nargin == 0
        sigmainit = 0.9
 %       sigmainit = [sigmainit; sigmainit];
    else
        sigmainit = varargin{1};
    end
tic
display('*****  Estimating  *****');

%for ces = 0:1
ces = 1

RandStream.setGlobalStream( RandStream('mt19937ar','Seed', 1));
%load simmarket;

%% ESTIMATION OF DEMAND

if ces == 1
    load simdataces;
    demand = MixedLogitDemand(nevo);
    demand.settings.ces = true;
else
    load simdata;
    demand = MixedLogitDemand(nevo);
end
demand.var.exog = 'd';
demand.var.endog = [ ];
% demand.var.panel = 'brn';
demand.var.market = 'marketid';
demand.var.price = 'p';
demand.settings.paneltype = 'none';
demand.var.quantity = 'quantity';
demand.var.marketsize = 'ms';
demand.var.instruments = 'inst1 inst2 inst3 inst4 inst5 inst6';
%demand.instruments = 'inst1 inst2 inst3';
sigmainit
demand.rc_sigma =  sigmainit;      % Gives initial guess for sigma's  
demand.var.nonlinear = 'd';

demand.settings.nind = 300;
demand.settings.drawmethod = 'hypercube';


demand.init();

res1= demand.gls();
res1= demand.gmm()

results = demand.estimate();
if false
    demand.optimalIV = true;
    results = demand.estimate();
end
results.thetaorg = betatrue;
disp(results);
toc
