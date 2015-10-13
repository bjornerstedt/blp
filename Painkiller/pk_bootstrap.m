%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RC bootstrap with cost and demand calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tic
% Primary Config parameters

RC = true    % Set RC to choose between RC and NL
ces = true   % Set ces to choose between ces and unit demand
select = true  % Choose between multiple RC starting points or selection  

% Other parameters

estimate = true % otherwise use saved estimate
display = true;
saveXLS = true;
bootreps = 1
conduct = [0, 0.75];

input.repetitions = 1;
input.save = true;
input.newinstruments = true;
input.withtime = false;
input.drawmethod = 'quadrature';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if RC
    if ces
        if input.withtime
        fn = 'pk RC CES time boot';
        selectRun = [1, 3]; % CES demand
        else
        fn = 'pk RC CES boot quad';
        selectRun = [1, 1]; % CES demand
        end
        dfn = 'demandCES'; % Saved simulation 3
    else
        selectRun = [3, 1]; % Unit demand
        fn = 'pk RC Unit boot quad';
        dfn = 'demandUnit'; % Saved simulation 7
    end
    if estimate
        nonlinear = cell(4,1);
        nonlinear(:) = {'paracetamol fizzy branded constant'};
        nind  = [500; 500; 500; 500];
        cesdemand = {ces; ces; ces; ces};
        quaddraws = [10; 11; 12; 13];

        testdata = table2struct(table(nonlinear, cesdemand, nind, quaddraws));
  %      testdata(2, :) = [];
        if select
            input.nocons = true;
            input.packDemand = false;
            [job,diary] = batchrun(@pkRC, input, testdata, input.repetitions, ...
                'Parallel', false, 'select', selectRun, ...
                'randomstream', 99);
            demand = job{1}.demand;
         %   save(dfn, 'demand') % Saved selected simulation
        else
            input.nocons = false;
            input.fn = 'pk RC Quadrature comp';
            input.packDemand = true;
            [job,diary] = batchrun(@pkRC, input, testdata, input.repetitions, ...
                'process', @pkRC_output, 'Parallel', false, 'randomstream', 99);
            toc
            break % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        end
    else
        load(dfn) % Saved simulation 3
    end
else
    if ces
        fn = 'pk NL CES boot';
    else
        fn = 'pk NL Unit boot';
    end
    pkNL % Estimate NL demand, uses ces parameter
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newdemand = copy(demand);

selection = (demand.T.year==2008 & demand.T.month==12);

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
bootpricech = zeros(6, bootreps); % length(unique(market.firm))
bootcosts = zeros(6, bootreps); % length(unique(market.firm))
disp('Parametric Bootstrap of Market Results')
names = cell(bootreps,1);
meanpc = [];
tic
load alldataUNIT12;
% alldata = table;

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
            newdemand.beta(1:size(newdemand.results.estimate, 1)); % HACK
    end
end

newdemand.init();

vars = strsplit('product year month ldosage lpacksize form1 substance1 substance2 branded');
proddata = demand.T(:, vars);
proddata = sortrows(proddata, {'product', 'year'},{'ascend','descend'});
diff = [0; proddata.product(1:end-1)] == proddata.product;
proddata(diff,:) = [];
proddata.product = [];
proddata2 = proddata;
proddata.delta = delta;
%proddata(proddata.year ~= 2008,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% proddata2(proddata2.year ~= 2008,:) = [];

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
% Display
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
%    xlswrite(filename, {'Second stage estimates - without and with collusion'}, 'Estimate', sprintf('A%d', pos-3));
    et.setColwise(outcol{cc});
    et.write(demandreg.results.estimate, 'Heading', 'Demand Estimate', 'RowNamesHeading', 'Estimate');
    et.write(costreg.results.estimate, 'Heading', 'Cost Estimate', 'RowNamesHeading', 'Estimate');
    et.write(market1.summary( 'brand'), 'Heading', 'Original market', 'RowNamesHeading', 'brand');
    et.write(market1.summary( 'substance'), 'RowNamesHeading', 'substance');
    et.write(marketNC.summary( 'brand'), 'Heading', 'New costs market', 'RowNamesHeading', 'brand');
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
if bootreps > 1
    save(fn, 'alldata');
    pk_bootstrap_output
end
toc