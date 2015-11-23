%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP MERGER SIMULATION WITH PAINKILLER DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametric bootstrap merger simulation
clear;
estimateDemand = false
nonlinp = false
optimalIV = true
quadrature = false
merge = true
newinstruments = false
drawBeta = false            % Study effect of individual draws in RC
bootreps = 1

display = true;
saveXLS = false;
ces = 1;

stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

%load painkillers9511new;
load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');
pk.constant = ones(size(pk,1),1);

if ces == 1
    fn = 'paracetamolCESML';
    pk.Xtablets = pk.Xtablets*10e-7;
    demand = CesMixedLogitDemand(pk);
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
else
    fn = 'paracetamolUnitML';
    demand = MixedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
end

demand.var.quantity = 'Xtablets';
demand.var.market = 'time';
demand.var.panel = 'product';
demand.var.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];

if newinstruments
demand.var.instruments = ['i1_con i2_con i1_ldosage i2_ldosage i1_lpacksize '...
    'i2_lpacksize i1_form2 i2_form2 i1_substance2 i2_substance2 i1_substance3 i2_substance3'];
else
demand.var.instruments = 'num numg numf numfg numhg numfgh';
end

demand.settings.drawmethod = 'hypercube';
demand.settings.marketdraws = true;
demand.settings.nind = 1000;
demand.settings.paneltype = 'lsdv';
demand.settings.quaddraws = 10;
% demand.settings.parallel = true;
% demand.settings.guessdelta = true;
%demand.settings.fptolerance1 = 1e-8;

if nonlinp
%    demand.var.nonlinear = 'paracetamol ibuprofen asa Ptablets';
    demand.var.nonlinear = 'paracetamol Ptablets'; % -6.2892 6.1536
   % demand.var.nonlinear = 'Ptablets paracetamol';  % -8.3531 15.764
%demand.analyticJacobian = false; % Share jacobian in equilibrium finding
%demand.rc_sigma = 5;
else
    demand.var.nonlinear = 'paracetamol ibuprofen asa';
  %  demand.var.nonlinear = 'paracetamol branded';
    demand.var.nonlinear = 'paracetamol fizzy branded constant';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optionally use saved first stage estimation
if estimateDemand
    demand.init();
    if optimalIV
        demand.estimate();
        demand.settings.optimalIV = true;
    end
    results = demand.estimate();
    
    beta = demand.beta;
    xi = demand.xi;
    rc_sigma = demand.rc_sigma;
    edelta = demand.edelta;
    varcovar = demand.varcovar;
    alpha = demand.alpha;
        save estimate_MLboot beta alpha rc_sigma xi results varcovar;
else
        load estimate_MLboot;
    demand.beta = beta;
    demand.alpha = alpha;
    demand.rc_sigma = rc_sigma;
   % demand.xi = xi;
    demand.varcovar = varcovar;
    demand.results.estimate = results;
    results
end
    demand.init();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Standard Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newdemand = copy(demand);
newdemand.init(); 

selection = (pk.year==2008 & pk.month==12);

newdemand.initSimulation(selection);

market = Market(newdemand);
market.var.firm = 'firm';
market.init();
market.findCosts();

market2 = copy(market);

market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
% market2.firm = 2; % Monopoly
market2.p0 = market.p;
market2.init(); % Initialize ownership matrix again
market2.equilibrium();

[nlMergerResult, mpc] = market.compare(market2.p);
disp 'Merger result on estimated values'
nlMergerResult
mpc
market2.equilibrium();
[nlMergerResult, mpc] = market.compare(market2.p);
disp 'Merger result on estimated values'
nlMergerResult
mpc

newdemand2 = copy(demand);
newdemand2.init(); 

selection = (pk.year==2008 & pk.month==12);

newdemand2.initSimulation(selection);

market = Market(newdemand);
market.var.firm = 'firm';
market.init();
market.findCosts();

market2 = copy(market);

market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
% market2.firm = 2; % Monopoly
market2.p0 = market.p;
market2.init(); % Initialize ownership matrix again
market2.equilibrium();

[nlMergerResult, mpc] = market.compare(market2.p);
disp 'Merger result on estimated values'
nlMergerResult
mpc

