%% Test 1: Paracetamol 

clear;
tic
allTests = false 
test = 1
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
demand = RCDemand(pk);
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
demand.settings.paneltype = 'lsdv';
demand.settings.marketdraws = true;

if test == 1 || allTests
dtest = copy(demand);

dtest.var.nonlinear = 'paracetamol';  % prices added if nonlinearPrice.
dtest.rc_sigma = [ -0.7648];      % Initial guess with wrong result  
dtest.rc_sigma = [ 5];      % Gives initial guess for sigma's  
dtest.init();

results = dtest.estimate();
assert(abs(results{'Ptablets','Coef'} + 6.601056128324467)<1e-4)
assert(abs(results{'rc_paracetamol','Coef'} - 5.094421916584797)<1e-2)
if ~strcmpi(demand.settings.paneltype, 'lsdv')
    assert(abs(results{'Ptablets','Std_err'} - 0.504236078715599)<1e-2)
    assert(abs(results{'rc_paracetamol','Std_err'} - 2.482284107760600)<1e-2)
end
display '****************** Test 1 passed **********************'

% 1b. Paracetamol with Optimal IV
dtest.settings.optimalIV = true;

results = dtest.estimate();

if strcmpi(demand.settings.paneltype, 'lsdv')
    assert(abs(results{'Ptablets','Coef'} + 6.276763124690481)<10e-4)
    assert(abs(results{'rc_paracetamol', 'Coef'} - 1.617406824151684)<1e-3)
else
    assert(abs(results{'Ptablets','Coef'} + 6.403504879090317)<10e-4)
    assert(abs(results{'rc_paracetamol', 'Coef'} - 3.081800213755897)<1e-3)
    assert(abs(results{'Ptablets','Std_err'} - 0.543038721313100)<10e-4)
    assert(abs(results{'rc_paracetamol','Std_err'} - 1.010761984537725)<1e-4)
end
display '****************** Test 1b passed **********************'
end

% Test 2: Constant 

if test == 2 || allTests
dtest = copy(demand);

dtest.var.nonlinear = 'constant';  % prices added if nonlinearPrice.
%dtest.rc_sigma = [ 5;.1];      % Gives initial guess for sigma's  
dtest.init();

results = dtest.estimate();

%assert(abs(results{'rc_constant', 'Coef'} + 0.4833)<1e-2)
%assert(abs(results{'rc_constant', 'Std_err'} - 5.6928)<1e-1)
%Compared to painkillBLP with exact share calc rather than from logshares
assert(abs(results{'rc_constant', 'Coef'} - 0.482326010258125)<10e-3)
assert(abs(results{'rc_constant', 'Std_err'} - 5.722081662730967)<10e-2)
display '****************** Test 2 passed **********************'
end
% Test 3: packsize 

if test == 3 || allTests

dtest = copy(demand);

dtest.var.nonlinear = 'packsize';  % prices added if nonlinearPrice.
dtest.rc_sigma = [ .1];      % Gives initial guess for sigma's  
dtest.init();

results = dtest.estimate();
assert(abs(results{'rc_packsize','Coef'} - 0.0017)<1e-4)
assert(abs(results{'rc_packsize','Std_err'} - 0.2161)<1e-3)
display '****************** Test 3 passed **********************'
end
% Test 4: Paracetamol and packsize

if test == 4 || allTests

dtest = copy(demand);

dtest.var.nonlinear = 'paracetamol packsize';  % prices added if nonlinearPrice.
dtest.rc_sigma = [ 5;.1];      % Gives initial guess for sigma's  
dtest.init();

results = dtest.estimate();

assert(abs(results{'rc_paracetamol','Coef'} - 5.08584575286999)<1e-2)
assert(abs(results{'rc_packsize','Coef'} - 0.000205876247160164)<1e-4)
display '****************** Test 4 passed **********************'
end
% Test 5: Paracetamol and price

if test == 5 || allTests

dtest = copy(demand);

dtest.var.nonlinear = 'paracetamol Ptablets';  % prices added if nonlinearPrice.
dtest.rc_sigma = [ .1;5];      % Gives initial guess for sigma's  
%dtest.exponentialFPiteration = true;
dtest.init();
results = dtest.estimate();
% Results from painkillBLP with Xrandom=[paracetamol Ptablets ]
% NOT the same as in resultsP.xls
assert(abs(results{'rc_paracetamol', 'Coef'} - 16.6025)<1e-3)
assert(abs(results{'rc_Ptablets', 'Coef'} - 7.3369)<1e-3)
assert(abs(results{'rc_paracetamol', 'Std_err'} - 3.4042)<1e-1)
assert(abs(results{'rc_Ptablets', 'Std_err'} - 0.8277)<1e-2)
display '****************** Test 5 passed **********************'
if false
    dtest.settings.optimalIV = true;
    results = dtest.estimate();
end
end
toc
