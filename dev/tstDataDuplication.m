clear
getpainkiller

pk.paracetamol = +(pk.substance =='Paracetamol');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
d1 = RCDemand(pk);
d2 = RCDemand(pk);

% Define demand
d1.var.quantity = 'Xtablets';
d1.var.marketsize = 'BL_Unit';
d1.var.market = 'time';
d1.var.panel = 'product';
d1.var.exog = ['marketing1 sw sm  month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
d1.var.price = 'Ptablets'; 
d1.var.instruments = 'num numg numf numfg numhg numfgh';
d1.var.nonlinear = 'paracetamol'

d1.estimate()

d2 = copy(d1);
d1.data.ibuprofen = +(pk.substance =='Ibuprofen');
