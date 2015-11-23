%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMM Estimate of both demand and supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

ces = 0
mixedLogit = 0

load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mixedLogit == 0
    if ces == 1
        demand = CesNestedLogitDemand(pk);
    else
        demand = NestedLogitDemand(pk);
    end
    demand.nests = 'form substance';
    demand.exog = ['marketing1 sw sm  month2 month3 month4 month5 month6 '...
        'month7 month8 month9 month10 month11 month12'];
else
    if ces == 1
        demand = CesMixedLogitDemand(pk);
    else
        demand = MixedLogitDemand(pk);
    end
    demand.drawmethod = 'hypercube';
    demand.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
        'month7 month8 month9 month10 month11 month12'];
    demand.nonlinear = 'paracetamol';

end
if ces == 1
    pk.Xtablets2 = pk.Xtablets*10e-7;
    demand.marketsize = 'BL_CES';
    demand.price = 'Ptablets'; 
    demand.quantity = 'Xtablets2';
else
    demand.marketsize = 'BL_Unit';
    demand.price = 'Ptablets_Real'; 
    demand.quantity = 'Xtablets';
end
demand.market = 'date';
demand.panel = 'product';
demand.instruments = 'num numg numf numfg numhg numfgh';

demand.settings.paneltype = 'lsdv';
demand.nocons = true;
demand.init(); 
if mixedLogit == 0
    demand.gls()
    delta = demand.beta(size(demand.results,1)+1:end);
else
    demand.estimate();
    delta = demand.beta(size(demand.results,1):end); % ONE random coefficient!
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cost estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mixedLogit == 0
    demand.initSimulation();
else 
    demand.init();
end
market0 = GmmMarket(demand);
market0.market = 'date';
market0.panel = 'product';
market0.exog = ['time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
market0.settings.paneltype = 'lsdv';
market0.nocons = true;
market0.firmvar = 'firm';
market0.costfunction = 'loglinear';
market0.estimate();

[betaD betaM] = market0.minimize();
demand.results
market0.results
[betaD(1:10) betaM(1:10)]
