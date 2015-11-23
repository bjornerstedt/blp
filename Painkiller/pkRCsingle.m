%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST Mixed Logit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testMLpk.m tests mixed logit demand on painkiller dataset with
% * CES and unit demand
% Only unit demand works because of ill conditioned matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Load data and prepare demand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
tic
% stream1 = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(stream1);

optimalIV = true
startingPoints = 1
merge = true
parallel = false
ces = 1

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
    demand = MixedLogitDemand(pk);
    demand.settings.ces = true;
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
else
    demand = MixedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
end

demand.var.nonlinear = 'paracetamol'; 
% demand.var.nonlinear = 'paracetamol constant';

demand.var.quantity = 'Xtablets';
demand.var.market = 'time';
demand.var.panel = 'product';
demand.var.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.settings.paneltype = 'lsdv';
demand.var.instruments = 'num numg numf numfg numhg numfgh';

demand.settings.nind = 300;

demand.config.parallel=parallel; % Replicability (below) requires serial execution
% demand.settings.fptolerance1 = 1e-12; % use lower tolerance for first FP iterations
% demand.settings.fptolerance2 = 1e-12; % use maximum tolerance for last iterations
demand.settings.drawmethod = 'hypercube';
demand.settings.marketdraws = true;
demand.init(); 
    if optimalIV
        demand.estimate();
        demand.settings.optimalIV = true;
    end
% demand.config.startingPoints = startingPoints;
demand.estimate();
% demand.estimate('GradObj','off');
rc_sigmas0 = cell(0);
rc_sigmas = cell(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Mixed Logit Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~merge
    break;
end
%selection = (pk.year==2008 & pk.month==12);
selection = (pk.year==2008 );

demand.initSimulation(selection);

market = Market(demand);
market.var.firm = 'firm';
market.init();
market.findCosts();

market2 = copy(market);

market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
%market2.firm = 2*ones(size(market2.firm)); % Monopoly
market2.p0 = market.p;
market2.init(); % Initialize ownership matrix again

market2.equilibrium();
%market2.fixedPoint(1000)
nlMergerResult = market.compare(market2.p);
nlMergerResult
toc
