%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimentation file for ML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Load data and prepare demand
clear;
%parpool;
tic
estimateDemand = true
nonlinp = false
optimalIV = false
quadrature = false
merge = true

stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

load painkillers9511new;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');

% Define demand
demand = MixedLogitDemand(pk);
demand.quantity = 'Xtablets';
demand.marketsize = 'BL';
demand.market = 'time';
demand.panel = 'product';
demand.exog = ['marketing1 sw sm  month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.price = 'Ptablets'; 
demand.instruments = 'num numg numf numfg numhg numfgh';
demand.drawmethod = 'halton';
demand.nind = 300;
demand.settings.paneltype = 'lsdv';
demand.marketdraws = false;
demand.quaddraws = 10;
demand.exponentialFPiteration = true;
demand.parallel = true;
demand.guessdelta = true;
%demand.fptolerance1 = 1e-8;

if nonlinp
%    demand.nonlinear = 'paracetamol ibuprofen asa Ptablets';
    demand.nonlinear = 'paracetamol Ptablets'; % -6.2892 6.1536
   % demand.nonlinear = 'Ptablets paracetamol';  % -8.3531 15.764
%demand.analyticJacobian = false; % Share jacobian in equilibrium finding
%demand.rc_sigma = 5;
else
%    demand.nonlinear = 'paracetamol ibuprofen asa';
 %   demand.nonlinear = 'paracetamol branded';
    demand.nonlinear = 'paracetamol';
  %  demand.nonlineartriangular = 'paracetamol';
    
end

if false
demand.drawmethod = 'quadrature';

demand.quadrature = true;
demand.quaddraws = 5;
demand.nind = 300;
demand.init();
demand.estimate();
demand.draws = [];
demand.drawmethod = 'halton';
end
demand.quaddraws = 15;

demand.init();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optionally use saved first stage estimation
if estimateDemand
    if optimalIV
        demand.estimate();
        demand.optimalIV = true;
    end
    results = demand.estimate();
    
    beta = demand.beta;
    xi = demand.xi;
    rc_sigma = demand.rc_sigma;
    edelta = demand.edelta;
    if nonlinp
        save estimate_NLp beta rc_sigma xi results;
    else
        save estimate_NLsubst beta rc_sigma xi results;
    end
else
    if nonlinp
        load estimate_NLp ;
    else
        load estimate_NLsubst beta rc_sigma xi results;
    end
    demand.beta = beta;
    demand.rc_sigma = rc_sigma;
    demand.xi = xi;
    results
end
toc
%demand.rc_sigma = [.3;.1].*rc_sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~merge
    break
end
%sel = (pk.year==2008 & pk.month==12);
sel = (pk.year==2008);

demand.init([], sel);
market = Market(demand);
market.firmvar = 'firm';
market.init();
market.findCosts();

market2 = copy(market);

market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; % Buyer(2) seller(1)
% market2.firm = 2; % Monopoly
market2.p0 = market.p;
market2.init(); % Initialize ownership matrix again

market2.equilibrium();
%market2.fixedPoint(1000)
[mergerresult,apc] = market.compare(market2.p);
mergerresult
apc


toc
% if ~quadrature
% mergerresultFE = mergerresult;
% resultsFE = results
% save temp resultsFE mergerresultFE
% else
%     load temp
%     results
%     resultsFE
%     mergerresult
%     mergerresultFE
% end
