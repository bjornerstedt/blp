%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP MERGER SIMULATION WITH PAINKILLER DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_painkiller.m replicates merger results of 
% C:\Users\n13017\Documents\Arbete\Painkiller\matlab_comparison.do

%% 0. Load data and prepare demand
clear;

load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Nested Logit Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define nested logit demand

    demand = NestedLogitDemand(pk);
    demand.marketsize = 'BL_Unit';
    demand.price = 'Ptablets_Real'; 

demand.nests = 'form substance';
demand.quantity = 'Xtablets';
demand.market = 'date';
demand.panel = 'product';
demand.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.instruments = 'num numg numf numfg numhg numfgh';

demand.settings.paneltype = 'lsdv';

demand.init(); 
demand.gls()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Nested Logit Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selection = (pk.year==2008 & pk.month==12);
demand.initSimulation(selection);

market = Market(demand);
market.firmvar = 'firm';
market.init();
market.findCosts();

market2 = copy(market);
market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
% market2.firm = 2; % Monopoly
market2.p0 = market.p;
market2.init(); % Initialize ownership matrix again

market2.equilibrium();
%market2.fixedPoint(1000)
mergerResult = market.compare(market2.p);
mergerResult

assert(abs(mergerResult{1,'Price2'} - 0.533598  )<10e-3)
assert(abs(mergerResult{1,'PriceCh'} - 0.1059824  )<10e-3)
display '************** Mergersim Unit Test passed ******************'

A = demand.G;
S = demand.s;
P = market.p;
D = demand.shareJacobian(market.p);
    
grelas = A*diag(P) * D'* A' * diag( 1 /(A*S) ) 
elas = diag(P)  * D' *diag( 1 / S );

