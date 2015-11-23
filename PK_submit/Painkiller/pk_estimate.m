%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RC bootstrap with cost and demand calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% Primary Config parameters
RC = true    % Set RC to choose between RC and NL
ces = true   % Set ces to choose between ces and unit demand

% Other parameters
show = true;
saveXLS = true;
bootreps = 1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stream1 = RandStream('mt19937ar','Seed',99);
RandStream.setGlobalStream(stream1);
load painkillers;
fnames = {'estimate_unit','estimate_CES'};
fn = fnames{ces+1};
if RC
    disp '**************************************************************'
    disp '******************  Mixed Logit Demand  **********************'
    disp '**************************************************************'
    disp ' '
    demand = MixedLogitDemand(pk(pk.year<=2008, :));
    if ces
        disp 'CES Demand';
        demand.settings.ces = true;
        demand.var.marketsize = 'BL_CES';
        demand.var.price = 'Ptablets';
        demand.var.quantity = 'Xtablets1';
    else
        disp 'Unit Demand';
        demand.var.marketsize = 'BL_Unit';
        demand.var.price = 'Ptablets_Real';
        demand.var.quantity = 'Xtablets';
    end
    
    demand.var.nonlinear = 'paracetamol fizzy branded constant';
    
    demand.var.market = 'time';
    demand.var.panel = 'product';
    
    demand.var.exog = ['marketing sw sm month2 month3 month4 month5 month6 '...
        'month7 month8 month9 month10 month11 month12'];
    demand.var.instruments = ['i1_con i2_con i1_ldosage i2_ldosage i1_lpacksize '...
        'i2_lpacksize i1_form2 i2_form2 i1_substance2 i2_substance2 '...
        'i1_substance3 i2_substance3'];
    
    demand.settings.drawmethod = 'quadrature';
    demand.settings.paneltype = 'lsdv';
    demand.settings.quaddraws = 10;
    demand.settings.nocons = true;
    demand.init();
    
    demand.estimate();
    optimalIV = true;
    if optimalIV
        disp('****** Optimal IV estimation ******');
        demand.settings.optimalIV = true;
        demand.estimate();
    end
else
    if ces
        fn = 'pk NL CES boot';
    else
        fn = 'pk NL Unit boot';
    end
    disp '**************************************************************'
    disp '*****************  Nested Logit Demand  **********************'
    disp '**************************************************************'
    disp ' '
    demand = NestedLogitDemand(pk(pk.year<=2008, :));
    demand.settings.ces = ces;
    demand.var.nests = 'form substance';
    demand.var.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
        'month7 month8 month9 month10 month11 month12'];
    demand.config.quietly = true;
    if ces
        disp 'CES Demand'
        demand.var.marketsize = 'BL_CES';
        demand.var.price = 'Ptablets';
        demand.var.quantity = 'Xtablets1';
    else
        disp 'Unit Demand'
        demand.var.marketsize = 'BL_Unit';
        demand.var.price = 'Ptablets_Real';
        demand.var.quantity = 'Xtablets';
    end
    demand.var.market = 'date';
    demand.var.panel = 'product';
    demand.var.instruments = 'num numg numf numfg numhg numfgh';
    
    demand.settings.paneltype = 'lsdv';
    demand.settings.nocons = true;
    demand.init();
    demand.estimate();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selection = (demand.T.year==2008 & demand.T.month==12);
newdemand = copy(demand);
newdemand.initSimulation(selection);

market = Market(newdemand);
market.var.firm = 'firm';
market.init();
market.findCosts();

disp 'Elasticities in initial market';
elsubst = newdemand.groupElasticities(market.p, 'substance')
elbrand = newdemand.groupElasticities(market.p, 'brand')
if RC
    [elas, E] = newdemand.elasticities(market.p, 'substance');
else
    [elas, E] = newdemand.elasticities(market.p);
end
if saveXLS
    et = ExcelTable([fn,'.xlsx']);
    et.setSheet('Estimate');
    et.write(demand.results.settings, 'Heading', 'Settings', 'RowNamesHeading', 'Settings');
    et.write(demand.results.var, 'Heading', ' ', 'RowNamesHeading', 'Variable',  'below', false);
    et.write(demand.results.estimate, 'Heading', 'Demand Estimate', 'RowNamesHeading', 'Estimate');
    et.write(elas, 'Heading', 'Product Elasticities', 'RowNamesHeading', 'Elasticity');
    et.write(elbrand, 'Heading', 'Group Elasticities', 'RowNamesHeading', 'brand');
    et.write(elsubst, 'Heading', 'Group Elasticities', 'RowNamesHeading', 'substance');
    xlswrite([fn,'.xlsx'], E, 'Elasticities'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(newdemand); 

delta = demand.results.betadummies;

disp 'Demand estimate'
disp(demand.results.estimate);

selection = (demand.T.year==2008 & demand.T.month==12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RC
    sel = [{demand.getPriceName()}, demand.vars2];
else
    sel = [{demand.getPriceName()}, demand.getLogShareNames()]; 
end

varcovar = table2array( demand.results.params.varcovar(sel,sel));
varcovar = (varcovar + varcovar')/2;
bootpars = mvnrnd(...
    table2array( demand.results.estimate(sel,'Coef')),...
    varcovar, bootreps);
bootdraws = array2table(bootpars);
bootdraws.Properties.VariableNames = sel;
bootresults = zeros(length(demand.vars), bootreps);
bootpricech = zeros(6, bootreps); 
bootcosts = zeros(6, bootreps); 
disp('Parametric Bootstrap of Market Results')
names = cell(bootreps,1);
meanpc = [];

conduct = [0, 0.75];
alldata = table;
parfor_progress(bootreps);
for r = 1:bootreps
for cc = 1:length(conduct)    

newdemand = copy(demand);
if r > 1 % First draw is the estimated value
    if RC
        newdemand.beta(1) =  bootpars(r, 1);
        newdemand.rc_sigma = bootpars(r, 2:end);
    else
        newdemand.beta(1:1+length(newdemand.nestlist)) = bootpars(r, :)';
        newdemand.results.estimate{:,1} = ...
            newdemand.beta(1:size(newdemand.results.estimate, 1)); 
    end
end
newdemand.init();

vars=strsplit('product year month ldosage lpacksize form1 substance1 substance2 branded');
proddata = demand.T(:, vars);
proddata = sortrows(proddata, {'product', 'year'},{'ascend','descend'});
diff = [0; proddata.product(1:end-1)] == proddata.product;
proddata(diff,:) = [];
proddata.product = [];
proddata2 = proddata;
proddata.delta = delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand delta estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demandreg = Estimate(proddata);
demandreg.var.exog = 'ldosage lpacksize form1 substance1 substance2 branded';
demandreg.settings.paneltype  = 'none';
demandreg.var.depvar = 'delta';
demandreg.init();
demandest = demandreg.estimate();
disp 'Delta estimation'
disp(demandest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cost estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demand0 = copy(newdemand);
demand0.initSimulation();

market0 = Market(demand0);
market0.settings.conduct = conduct(cc);

market0.var.market = 'date';
market0.var.panel = 'product';
market0.var.exog = ['time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
market0.settings.paneltype = 'lsdv';
market0.settings.nocons = true;
market0.var.firm = 'firm';
market0.costfunction = 'loglinear';
market0.estimateCosts();

costs = market0.results.betadummies;
proddata2.costs = costs;

costreg = Estimate(proddata2);
costreg.var.exog = 'ldosage lpacksize form1 substance1 substance2';
costreg.settings.paneltype  = 'none';
costreg.var.depvar = 'costs';
costreg.init();
costest = costreg.estimate();
disp 'Cost estimate'
disp(costest);

elim = ( demand.T.elim == 0);
elim = elim(selection);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merger/market Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newdemand.initSimulation(selection);

delta_new =  newdemand.d + ...
    (log(newdemand.T.packsize_new) - newdemand.T.lpacksize) .* ...
    demandest{'lpacksize', 'Coef'};
costs_new = market0.c(selection) .* ...
    (newdemand.T.packsize_new ./ newdemand.T.packsize) .^ ...
    costest{'lpacksize', 'Coef'};

market1 = Market(newdemand);
market1.var.firm = 'firm';
market1.settings.conduct = conduct(cc);
market1.init();
market1.c = market0.c(selection);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Computing market equilibrium with new costs'

marketNC = copy(market1);
marketNC.p0 = market1.p;
marketNC.init(); % Initialize ownership matrix again
marketNC.c = costs_new;
marketNC.equilibrium();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Computing market equilibrium with merger'

market2 = copy(market1);
market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
market2.p0 = market1.p;
market2.init(); % Initialize ownership matrix again
market2.equilibrium();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Computing market equilibrium with merger and new costs'

market2NC = copy(market2);
market2NC.p0 = market2.p;
market2NC.init(); % Initialize ownership matrix again
market2NC.c = costs_new;
market2NC.equilibrium();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Computing market equilibrium with merger, new costs and demand'

market2NCdem = copy(market2NC);
market2NCdem.D.d = delta_new;
market2NCdem.equilibrium();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put data in table

if conduct(cc) == 0
    markets = {market1, marketNC, market2, market2NC, market2NCdem};
    resultcat = {'Original' , 'NC' , 'Merge', 'MergeNC', 'MergeNCD'};
else
    markets = { market1, marketNC, market2, market2NC, market2NCdem};
    resultcat = { 'Coll' , 'NCcoll' , 'MergeColl', 'MergeNCcoll', 'MergeNCDcoll'};
end
nc = cell(size(market1.p));
const = ones(size(market1.p));
for i = 1:length(markets)
    sharetab = markets{i}.summary();
    sharetab.Result = nc;
    sharetab.Result(:) = resultcat(i);
    sharetab.Bootrep = const * r;
    sharetab.Firm2 = market1.T.firmsubst;
    sharetab.Substance = newdemand.T.substance;
    alldata = [alldata; sharetab];
end

if saveXLS && r == 1
    outcol = {'A', 'H'};
    et.setColwise(outcol{cc});
    et.write(demandreg.results.estimate, 'Heading', 'Demand Estimate', ...
        'RowNamesHeading', 'Estimate');
    et.write(costreg.results.estimate, 'Heading', 'Cost Estimate', ...
        'RowNamesHeading', 'Estimate');
    et.write(market1.summary( 'brand'), 'Heading', 'Original market', ...
        'RowNamesHeading', 'brand');
    et.write(market1.summary( 'substance'), 'RowNamesHeading', 'substance');
    et.write(marketNC.summary( 'brand'), 'Heading', 'New costs market', ...
        'RowNamesHeading', 'brand');
    et.write(marketNC.summary( 'substance'), 'RowNamesHeading', 'substance');
end
end
delhandles = {newdemand,demand0,demandreg,market0, market1, market2, costreg, ...
    marketNC, market2NC, market2NCdem};
for h = 1:length(delhandles)
    delete(delhandles{h}); 
end
parfor_progress;
end
parfor_progress(0);
alldata.Result = categorical(alldata.Result);
    
pk_bootstrap_output

