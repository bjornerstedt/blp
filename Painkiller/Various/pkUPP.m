%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP MERGER SIMULATION WITH PAINKILLER DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testUPP.m calculates UPP for multi-product market

%% 0. Load data and prepare demand
clear;

load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);

% subtotal = accumarray(pk.time, pk.PX1);
% subtotal =  subtotal(pk.time,:);
% subtotal = mean(subtotal);
% medincome = 649144.4; % time varying potential budget measure
% pk.BL = 2*subtotal * pk.GDPnom /medincome;	
% pk.quantity = pk.PX1 ./ pk.Ptablets;
%pk.marketing = pk.marketing*10e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Nested Logit Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define nested logit demand

    demand = NestedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
demand.var.nests = 'form substance';
demand.var.quantity = 'Xtablets';
demand.var.market = 'date';
demand.var.panel = 'product';
demand.var.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.var.instruments = 'num numg numf numfg numhg numfgh';

demand.settings.paneltype = 'lsdv';

demand.init(); 
demand.estimate()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Nested Logit Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selection = (pk.year==2008 & pk.month==12);
%selection = (pk.year==2008 );
firm = pk.firm(selection,:);
demand.initSimulation(selection);
market = Market(demand);
market.var.firm = 'firm';
market.init();
market.findCosts();

gsk = firm == 'GSK';
az = firm == 'AstraZeneca';
sj = demand.shareJacobian(demand.p);
div12 = -(gsk'*sj*az)/(gsk'*sj*gsk)
div21 = -(gsk'*sj*az)/(az'*sj*az)
e = .25;
% e = .1;
m = market.p - market.c;
meff = market.p - (1 - e)*market.c;
disp 'UPP 12'
upp12 = -market.c(gsk)*e - sj(gsk,gsk)\sj(gsk,az)*m(az)
mean(upp12)
-mean(market.c(gsk))*e + div12* mean(m(az)) 
disp 'UPP 21'
upp21 = -market.c(az)*e - sj(az,az)\sj(az,gsk)*m(gsk)
mean(upp21)
-mean(market.c(az))*e + div12* mean(m(gsk)) 
display 'UPP*12'
upp12 = -market.c(gsk)*e - sj(gsk,gsk)\sj(gsk,az)*meff(az)
mean(upp12)
upp21 = -market.c(az)*e - sj(az,az)\sj(az,gsk)*meff(gsk)
mean(upp21)

mc = mean(market.c)
