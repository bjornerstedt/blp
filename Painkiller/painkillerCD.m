%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RC and NL cost and demand calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost and demand delta change in nested logit

clear;
newinstruments = false
estimateDemand = false
optimalIV = false

for mixedLogit = 0:0
for ces = 1:1
    disp ''
    disp ''
    
load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');
pk.time = pk.date;
pk(pk.year>2008, :) = [];
pk.marketing = pk.marketing*10^-6;
pk.lpacksize = log(pk.packsize);
pk.ldosage = log(pk.dosage);
[~,~,pk.brandid]=unique(pk.brand);

pk.Xtablets2 = pk.Xtablets*10e-7;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ' ' 
disp ' ' 
if mixedLogit == 0
    disp '**************************************************************'
    disp '*****************  Nested Logit Demand  **********************'
    disp '**************************************************************'
    disp ' ' 
    if ces == 1
        demand = CesNestedLogitDemand(pk);
    else
        demand = NestedLogitDemand(pk);
    end
    demand.var.nests = 'form substance';
    demand.var.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
        'month7 month8 month9 month10 month11 month12'];
else
    disp '**************************************************************'
    disp '******************  Mixed Logit Demand  **********************'
    disp '**************************************************************'
    disp ' ' 
    if ces == 1
        demand = CesMixedLogitDemand(pk);
    else
        demand = MixedLogitDemand(pk);
    end
    demand.settings.drawmethod = 'hypercube';
    demand.var.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
        'month7 month8 month9 month10 month11 month12'];
    demand.var.nonlinear = 'paracetamol fizzy branded constant';
    demand.var.nonlinear = 'paracetamol';

end
if ces == 1
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
else
    demand.var.instruments = 'num numg numf numfg numhg numfgh';
end

demand.settings.paneltype = 'lsdv';
demand.settings.nocons = true;
demand.init(); 
if mixedLogit == 0
    results = demand.estimate();
    delta = demand.beta(size(results,1)+1:end);
else
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
    
    results = demand.estimate();
    delta = demand.beta(size(results,1)- length(demand.nonlinparams)+1:end ); 
end
disp 'Demand estimate'
disp(results);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars = strsplit('product year month ldosage lpacksize form2 substance2 substance3');
proddata = pk(:, vars);
proddata = sortrows(proddata, {'product', 'year'},{'ascend','descend'});
diff = [0; proddata.product(1:end-1)] == proddata.product;
proddata(diff,:) = [];
proddata.product = [];
proddata2 = proddata;
proddata.delta = delta;
%proddata(proddata.year ~= 2008,:) = [];

demandreg = Estimate(proddata);
demandreg.var.exog = 'ldosage lpacksize form2 substance2 substance3';
demandreg.settings.paneltype  = 'none';
demandreg.var.depvar = 'delta';
demandreg.init();
demandest = demandreg.estimate();
disp 'Delta estimation'
disp(demandest);

utilitychange = (log(pk.packsize_new)-pk.lpacksize).*demandest{'lpacksize', 'Coef'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cost estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mixedLogit == 0
    demand.initSimulation();
else 
    demand.init();
end
market0 = Market(demand);
market0.var.market = 'date';
market0.var.panel = 'product';
market0.var.exog = ['time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
market0.settings.paneltype = 'lsdv';
market0.settings.nocons = true;
market0.var.firm = 'firm';
market0.costfunction = 'loglinear';
market0.estimate();

costs = market0.betadummies;
proddata2.costs = costs;
% proddata2(proddata2.year ~= 2008,:) = [];

costreg = Estimate(proddata2);
costreg.var.exog = 'ldosage lpacksize form2 substance2 substance3';
costreg.settings.paneltype  = 'none';
costreg.var.depvar = 'costs';
costreg.init();
costest = costreg.estimate();
disp 'Cost estimate'
disp(costest);

costs_new = market0.c .* (pk.packsize_new ./ pk.packsize) .^ costest{'lpacksize', 'Coef'};

selection = (pk.year==2008 );
selection = (pk.year==2008 & pk.month == 12);
elim = ( pk.elim == 0);
elim = elim(selection);
%costChange = costChange(selection);
costs_newtest = costs_new(selection);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merger/market Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

demand.initSimulation(selection);

delta_new =  demand.d + utilitychange(selection,:);
market1 = Market(demand);
market1.var.firm = 'firm';
market1.init();
market1.c = market0.c(selection);

marketNC = copy(market1);
marketNC.p0 = market1.p;
marketNC.init(); % Initialize ownership matrix again
marketNC.c = costs_new(selection);

disp 'Computing market equilibrium with new costs'
marketNC.equilibrium();

market2 = copy(market1);
market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK'; 
% market2.firm = 2; % Monopoly
market2.p0 = market1.p;
market2.init(); % Initialize ownership matrix again

disp 'Computing market equilibrium with merger'
market2.equilibrium();

disp 'Computing market equilibrium with merger and new costs'
market2NC = copy(market2);
market2NC.p0 = market2.p;
market2NC.init(); % Initialize ownership matrix again
market2NC.c = costs_new(selection);
market2NC.equilibrium();

disp 'Computing market equilibrium with merger, new costs and new utility'
market2NCdem = copy(market2NC);
market2NCdem.D.d = delta_new;
market2NCdem.equilibrium();

ncResult = market1.compare(marketNC.p);
mergerResult = market1.compare(market2.p);
mergerResultNC = market1.compare(market2NC.p);
mergerResultNCdem = market1.compare(market2NCdem.p);

disp 'Price increase due to cost increase'
disp(ncResult);
disp 'Merger result without demand or cost changes'
disp(mergerResult);
disp 'Merger result with cost changes'
disp(mergerResultNC);
disp 'Merger result with demand or cost changes'
disp(mergerResultNCdem);

costChange = (market2NC.c - market1.c) ./ market1.c;

tableCols = {'Firm', 'Costs' 'Price1' 'Price2' 'PriceCh'};
tableCols = {'Firm', 'Costs' 'CostChange' };
tableColsNew = {'Count', 'Costs', 'Mean', 'Std'};

res = table(market1.firm, market1.c, ...
    costChange, 'VariableNames', tableCols);
res1 = varfun(@mean, res, 'GroupingVariables', 'Firm');
res2 = varfun(@std, res, 'GroupingVariables', 'Firm');
% Drop cols created by varfun
res2(:,'GroupCount') = []; 
res1(:,'Firm') = [];
res2(:,'Firm') = [];
res2(:,'std_Costs') = [];

disp 'Costs and percentage cost change by firm';
res1 = [res1 res2];
res1.Properties.VariableNames = tableColsNew;
disp(res1);

disp 'Costs and percentage cost change by firm for non-elim products';
res = res(elim,:);
res1 = varfun(@mean, res, 'GroupingVariables', 'Firm');
res2 = varfun(@std, res, 'GroupingVariables', 'Firm');
res2(:,'GroupCount') = [];
res1(:,'Firm') = [];
res2(:,'Firm') = [];
res2(:,'std_Costs') = [];
reselim = [res1 res2];
reselim.Properties.VariableNames = tableColsNew;
disp(reselim);
end
end