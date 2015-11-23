%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mixed Logit with Multiple starting points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test calculating exact deltas to see if there is a convergence problem in
% the contraction.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Load data and prepare demand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

optimalIV = true

ces = 0
draws = 2

load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Mixed Logit Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ces == 1
    pk.Xtablets = pk.Xtablets*10e-7;
    demandprefs = CesMixedLogitDemand(pk);
    demandprefs.marketsize = 'BL_CES';
    demandprefs.price = 'Ptablets'; 
else
    demandprefs = MixedLogitDemand(pk);
    demandprefs.marketsize = 'BL_Unit';
    demandprefs.price = 'Ptablets_Real'; 
end

demandprefs.quantity = 'Xtablets';
demandprefs.market = 'time';
demandprefs.panel = 'product';
demandprefs.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demandprefs.instruments = 'num numg numf numfg numhg numfgh';
demandprefs.settings.paneltype = 'lsdv';

demandprefs.settings.nind = 300;
demandprefs.nonlinearInstruments = false; 
demandprefs.exponentialFPiteration = false;

demandprefs.nonlinear = 'paracetamol'; 
demandprefs.settings.fptolerance1 = 1e-10; % use lower tolerance for first FP iterations
demandprefs.settings.fptolerance2 = 1e-14; % use maximum tolerance for last iterations
rc_sigmas0 = []; % Starting points
rc_sigmas1 = []; % Estimate
rc_sigmas2 = []; % Optimal IV
rc_sigmas3 = []; % Exact delta
for t=1:draws
    demand = copy(demandprefs);
    demand.init();
    demand.rc_sigma = randn(size(demand.rc_sigma));
    rc_sigmas0 = [rc_sigmas0, abs(demand.rc_sigma)];
    demand.estimate();
    rc_sigmas1 = [rc_sigmas1, abs(demand.rc_sigma)];
    
    demand.settings.optimalIV = true;
    demand.estimate();
    rc_sigmas2 = [rc_sigmas2, abs(demand.rc_sigma)];
    
    tic
    demand.exactDelta = true;
    demand.estimate();
    rc_sigmas3 = [rc_sigmas3, abs(demand.rc_sigma)];
    
    display(t);
    display 'Run time'
    toc
end
rc_sigmas0
rc_sigmas1
rc_sigmas2
rc_sigmas3