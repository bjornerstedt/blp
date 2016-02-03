function [results] = testMLsim(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP ESTIMATION WITH SIMULATED DATA SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETERMINATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ;
    if nargin == 0
        sigmainit = 0.8
 %       sigmainit = [sigmainit; sigmainit];
    else
        sigmainit = varargin{1};
    end
tic
display('*****  Estimating  *****');

for ces = 0:1
stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);
%load simmarket;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION OF DEMAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ces == 1
    load simdataces;
    demand = CesMixedLogitDemand(nevo);
else
    load simdata;
    demand = MixedLogitDemand(nevo);
end
demand.exog = 'f d';
demand.endog = [ ];
demand.panel = 'brn';
demand.market = 'cdid';
demand.price = 'p';
demand.settings.paneltype = 'none';
demand.quantity = 'quantity';
demand.marketsize = 'ms';
demand.instruments = 'inst1 inst2 inst3 inst4 inst5 inst6';
%demand.instruments = 'inst1 inst2 inst3';
sigmainit
demand.rc_sigma = [ sigmainit];      % Gives initial guess for sigma's  
demand.nonlinear = 'f';

demand.fpmaxit = 3000;
demand.nind = 300;
demand.fptolerance1 = 1e-12;
demand.analyticJacobian = true;

demand.init();

demand.gls();
demand.gmm()

results = demand.estimate();
if false
    demand.optimalIV = true;
    results = demand.estimate();
end
results.thetaorg = bethorg;
disp(results);

if ces  == 0
    assert(abs(results{'p','Coef'}  + 0.92728  )<10e-5)
    assert(abs(results{'rc_f','Coef'}  - 0.76155  )<10e-5)
    assert(abs(results{'rc_f','Std_err'}  - 0.17762  )<10e-5)
    display '************** Unit Demand Test passed ******************'
else
    assert(abs(results{'lP','Coef'}  + 0.79125  )<10e-5)
    assert(abs(results{'rc_f','Coef'}  - 0.6553  )<10e-5)
    assert(abs(results{'rc_f','Std_err'}  - 0.32634  )<10e-5)
    display '************** CES Demand Test passed ******************'
end
end