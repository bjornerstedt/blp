% Not yet a test, as it does not work

clear;
tic
allTests = true 
test = 1
optimalIV = true

stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

load painkillers9511new_orig;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');

% Define demand
demand = MixedLogitDemand(pk);
demand.var.quantity = 'Xtablets';
demand.var.marketsize = 'BL';
demand.var.market = 'time';
demand.var.panel = 'product';
demand.var.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.var.price = 'Ptablets'; 
demand.var.instruments = 'num numg numf numfg numhg numfgh';
demand.settings.drawmethod = 'halton2';
demand.settings.nind = 500;
demand.settings.marketdraws = true;

if test == 1 || allTests
dtest = copy(demand);

dtest.var.nonlinear = 'paracetamol';  % prices added if nonlinearPrice.
dtest.rc_sigma = [ -0.7648];      % Initial guess with wrong result  
dtest.rc_sigma = [ 5];      % Gives initial guess for sigma's  
dtest2 = copy(dtest);

dtest.settings.paneltype = 'fe';
dtest.init();

results = dtest.estimate();
assert(abs(results{'Ptablets','Coef'} + 6.601056128324467)<1e-4)
assert(abs(results{'rc_paracetamol','Coef'} - 5.094421916584797)<1e-2)
assert(abs(results{'Ptablets','Std_err'} - 0.504236078715599)<1e-2)
assert(abs(results{'rc_paracetamol','Std_err'} - 2.482284107760600)<1e-2)

if optimalIV
dtest.settings.optimalIV = true;
results = dtest.estimate();

assert(abs(results{'Ptablets','Coef'} + 6.403504879090317)<10e-4)
assert(abs(results{'rc_paracetamol', 'Coef'} - 3.081800213755897)<1e-3)
assert(abs(results{'Ptablets','Std_err'} - 0.543038721313100)<10e-4)
assert(abs(results{'rc_paracetamol','Std_err'} - 1.010761984537725)<1e-4)
end
display '****************** Test 1 FE passed **********************'
end
if test == 2 || allTests

dtest2.settings.paneltype = 'lsdv';
dtest2.init();

results = dtest2.estimate();
assert(abs(results{'Ptablets','Coef'} + 6.601056128324467)<1e-4)
assert(abs(results{'rc_paracetamol','Coef'} - 5.094421916584797)<1e-2)

assert(abs(results{'Ptablets','Std_err'} - 0.504236078715599)<1e-2)
assert(abs(results{'rc_paracetamol','Std_err'} - 2.482284107760600)<1e-2)

if optimalIV
dtest2.settings.optimalIV = true;
results = dtest2.estimate();

assert(abs(results{'Ptablets','Coef'} + 6.276763124690481)<10e-4)
assert(abs(results{'rc_paracetamol', 'Coef'} - 1.617406824151684)<1e-3)

%     assert(abs(results{'Ptablets','Coef'} + 6.403504879090317)<10e-4)
%     assert(abs(results{'rc_paracetamol', 'Coef'} - 3.081800213755897)<1e-3)
%     assert(abs(results{'Ptablets','Std_err'} - 0.543038721313100)<10e-4)
%     assert(abs(results{'rc_paracetamol','Std_err'} - 1.010761984537725)<1e-4)
end
display '****************** Test 1 LSDV passed **********************'
end
