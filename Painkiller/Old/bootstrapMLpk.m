%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP MERGER SIMULATION WITH PAINKILLER DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametric bootstrap merger simulation
clear;
estimateDemand = true
optimalIV = true
quadrature = false
merge = true
newinstruments = false
drawBeta = true            % Study effect of individual draws in RC
bootreps = 3

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
demand.settings.nind = 500;
demand.settings.paneltype = 'lsdv';
demand.settings.quaddraws = 10;
% demand.settings.parallel = true;
% demand.settings.guessdelta = true;
%demand.settings.fptolerance1 = 1e-8;

demand.var.nonlinear = 'paracetamol fizzy branded constant';
demand.var.nonlinear = 'paracetamol constant';

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
    newdemand = demand.pack();

    save 'test' 'newdemand';
    newdemand = copy(demand);
else
    load 'test'
    newdemand.init(pk); 
    demand = copy(newdemand);
    demand.results.estimate
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Standard Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selection = (pk.year==2008 & pk.month==12);

newdemand.init(pk, selection);

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

if ~merge
    break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Bootstrap Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sel = [{demand.getPriceName()}, demand.vars2];
varcovar = table2array( demand.results.params.varcovar(sel,sel));
varcovar = (varcovar + varcovar')/2;
bootpars = mvnrnd(...
    table2array( demand.results.estimate(sel,'Coef')),...
    varcovar, bootreps);

bootresults = zeros(length(demand.vars), bootreps);
bootpricech = zeros(6, bootreps); % length(unique(market.firm))
bootcosts = zeros(6, bootreps); % length(unique(market.firm))
disp('Parametric Bootstrap of Market Results')
names = cell(bootreps,1);
meanpc = [];
tic

parfor_progress(bootreps);
for r = 1:bootreps
% for r = 1:bootreps
% disp(r)
% Define mixed logit demand
newdemand = copy(demand);
if drawBeta
newdemand.beta(1) =  bootpars(r, 1);
newdemand.rc_sigma = bootpars(r, 2:end);
% newdemand.xi = [];
end
% newdemand.init(); 

selection = (pk.year==2008 & pk.month==12);

newdemand.initSimulation(selection);

nlMarket = Market(newdemand);
nlMarket.var.firm = 'firm';
nlMarket.init();
nlMarket.findCosts();

nlMarket2 = copy(nlMarket);

nlMarket2.firm(nlMarket2.firm == 'AstraZeneca' ) = 'GSK'; 
% market2.firm = 2; % Monopoly
nlMarket2.p0 = nlMarket.p;
nlMarket2.init(); % Initialize ownership matrix again
porig = nlMarket2.p0; 
for j = 1:5
    try
        nlMarket2.equilibrium('Display', 'off');
    catch err
        disp 'trying again'
        if j == 1
            nlMarket2.p0 = nlMarket.c;
        else
            nlMarket2.p0 = nlMarket.c + (porig - nlMarket.c)*rand;
        end
    end
end
%market2.fixedPoint(1000)
[nlMergerResult, mpc] = nlMarket.compare(nlMarket2.p);
if r == 1
    names{r} = nlMergerResult.Properties.RowNames;
end
bootpricech(:,r) = nlMergerResult{:,'PriceCh'};
bootcosts(:,r) = nlMergerResult{:,'Costs'};
meanpc = [meanpc mpc];
if display
    bootpars(r, :)
    disp(nlMergerResult)
end
parfor_progress;
end
parfor_progress(0);
cresultarray = [mean(bootcosts,2) std(bootcosts, 1, 2) min(bootcosts, [] ,2) max(bootcosts, [] ,2)];
presultarray = [mean(bootpricech,2) std(bootpricech, 1 ,2) min(bootpricech, [] ,2) max(bootpricech, [] ,2)];

headings = {'Coef', 'Std_err' , 'Min' , 'Max' };
cresults = array2table(cresultarray);
cresults.Properties.RowNames = names{1};
cresults.Properties.VariableNames = headings;

prresults = array2table(presultarray);
prresults.Properties.RowNames = names{1};
prresults.Properties.VariableNames =  headings;

disp('Starting beta');
bootpars'
disp('Costs');
cresults
disp('Price Change');
prresults

priceCh = array2table(bootpricech');
priceCh.Properties.VariableNames = names{1};
priceCh.meanpc = meanpc';
ci = sortrows(priceCh, 'meanpc');
ci = ci(ceil(0.05*bootreps), :);
disp 'Five percent one tail test'
ci

bootpricech = [bootpricech; meanpc]';

bootpricech = sortrows(bootpricech, size(bootpricech,2));
costs = array2table(bootcosts');
costs.Properties.VariableNames = names{1};
if saveXLS
xlswrite(fn, presultarray , 'Results', 'B3')
xlswrite(fn, cresultarray , 'Results', 'H3')
xlswrite(fn, names{1} , 'Results', 'A3')
xlswrite(fn, headings , 'Results', 'B2')
xlswrite(fn, names{1} , 'Results', 'G3')
xlswrite(fn, headings , 'Results', 'H2')

xlswrite(fn, {'Price Change'} , 'Results', 'A1')
xlswrite(fn, {'Costs' }, 'Results', 'G1')


xlswrite(fn, names{1}' , 'Price Change', 'A1')
xlswrite(fn, bootpricech , 'Price Change', 'A2')
xlswrite(fn, names{1}' , 'Costs', 'A1')
xlswrite(fn, bootcosts' , 'Costs', 'A2')
% xlswrite(fn, bootresults' , 'Estimate', 'A2')
% xlswrite(fn, demand.vars , 'Estimate', 'A1')
end
toc

%az = bootpricech(1,:)';
% hist(meanpc,50)
mean(priceCh{:,'meanpc'})