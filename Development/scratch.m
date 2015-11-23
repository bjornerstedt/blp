% Test condition of Z matrix when columns are dropped
% Very little gained

ces = true
optimalIV = true
newinstruments = true
withtime = false;
commonRC = 'constant' % RC in all treatments

load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');

if ces
    fn = 'paracetamolCESML'; 
    pk.Xtablets = pk.Xtablets*10e-7;
    demand = CesMixedLogitDemand(pk);
    demand.marketsize = 'BL_CES';
    demand.price = 'Ptablets'; 
else
    fn = 'paracetamolUnitML';
    demand = MixedLogitDemand(pk);
    demand.marketsize = 'BL_Unit';
    demand.price = 'Ptablets_Real'; 
end

demand.quantity = 'Xtablets';
demand.market = 'time';
demand.panel = 'product';
demand.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
if withtime
    demand.exog = ['time ' demand.exog];
end
if newinstruments
    demand.instruments = ['i1_con i2_con i1_ldosage i2_ldosage i1_lpacksize '...
        'i2_lpacksize i1_form2 i2_form2 i1_substance2 i2_substance2 '...
        'i1_substance3 i2_substance3'];
    % sum1_con sum2_con sum1_ldosage sum2_ldosage sum1_lpacksize sum2_lpacksize'
else
    demand.instruments = 'num numg numf numfg numhg numfgh';
end

%demand.drawmethod = 'halton';
demand.drawmethod = 'hypercube';
demand.marketdraws = false;
demand.nind = 300;
demand.settings.paneltype = 'lsdv';
demand.quaddraws = 10;
demand.fptolerance1 = 1e-12 ;
demand.fptolerance2 = 1e-14 ;

merger = Merger();
merger.selection = pk.year==2008 & pk.month==12;
merger.firmvar = 'firm';
merger.buyer = 'GSK';
merger.seller = 'AstraZeneca';

demand.nonlinear = ['paracetamol ',commonRC];

demand.init();

X = demand.X(:,2:end);
inst = table2array( pk(:,strsplit(strtrim(demand.instruments)) ));
Z=[X inst];
cond(Z'*Z)
c=zeros(11,11,11);
for i = 1:11
for j = 1:11
    if j==i continue; end
for k = 1:11
    if k==i continue; end
    a = sort([i,j,k]);
    inst1 = inst;
    inst1(:,a(3))= [];
    inst1(:,a(2))= [];
    inst1(:,a(1))= [];
    c(i, j, k)=cond([X inst1]);
end
end
end
min(min(min(c)))
max(max(max(c)))
cond([X inst])





break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results.estimate = demand.estimate();
results.merger = merger.merge(demand);
if optimalIV
    demand.optimalIV = true;
    results.estimate2 = demand.estimate();
    results.merger2 = merger.merge(demand);
end
delete(demand);
%delete(merger);



