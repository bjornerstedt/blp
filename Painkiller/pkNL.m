%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NL demand estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newinstruments = false
  
load painkillers9511main2new;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');
[~,~,pk.date] = unique(pk(:,{'year','month'}));
pk.time = pk.date + 419;
pk(pk.year>2008, :) = [];
pk.marketing = pk.marketing*10^-6;
pk.lpacksize = log(pk.packsize);
pk.ldosage = log(pk.dosage);
[~,~,pk.brandid]=unique(pk.brand);

pk.Xtablets2 = pk.Xtablets*10e-7;

pk.firm(pk.firm == 'Ellem') = 'Meda';
pk.firm(pk.firm == 'Recip') = 'Meda';
pk.firm(pk.firm == 'Pfizer') = 'McNeil';
pk.firmsubst = pk.firm;
pk.firmsubst(pk.brand == 'Ipren') = 'McNeil-Ibu';
pk.firmsubst(pk.brand == 'Alindrin') = 'Meda-Ibu';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);


disp '**************************************************************'
disp '*****************  Nested Logit Demand  **********************'
disp '**************************************************************'
disp ' ' 
demand = NestedLogitDemand(pk);
demand.settings.ces = ces;
demand.var.nests = 'form substance';
demand.var.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.config.quietly = true;

if ces
    disp 'CES Demand'
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
    demand.var.quantity = 'Xtablets2';
else
    disp 'Unit Demand'
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
    demand.var.quantity = 'Xtablets';
end
demand.var.market = 'date';
demand.var.panel = 'product';
if newinstruments
    demand.var.instruments = ['i1_con i2_con i1_ldosage i2_ldosage i1_lpacksize '...
        'i2_lpacksize i1_form2 i2_form2 i1_substance2 i2_substance2 '...
        'i1_substance3 i2_substance3'];
    demand.var.instruments = ['num numg numf numfg numhg numfgh' ...
        'i1_ldosage i1_lpacksize i2_ldosage i2_lpacksize'];
else
    demand.var.instruments = 'num numg numf numfg numhg numfgh';
end

demand.settings.paneltype = 'lsdv';
demand.settings.nocons = true;
demand.init(); 

results = demand.estimate();

