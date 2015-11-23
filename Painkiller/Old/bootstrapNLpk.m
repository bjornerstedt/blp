%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON-PARAMETRIC BOOTSTRAP NESTED LOGIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% 0. Load data and prepare demand
clear;
bootreps = 1;
display = false;
saveXLS = false
tic
for ces = 0:0

load 'painkillers9511main2new.mat';
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);

if ces == 1
    fn = 'cesresultsNL';
    pk.Xtablets = pk.Xtablets*10e-7;
    demand = CesNestedLogitDemand(pk);
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
else
    fn = 'unitresultsNL';
    demand = NestedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
end
demand.var.nests = 'form substance';
demand.var.quantity = 'Xtablets';
demand.var.market = 'date';
demand.var.panel = 'product';
demand.var.exog = ['marketing1 sw sm time month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.var.instruments = 'num numg numf numfg numhg numfgh';

demand.settings.paneltype = 'lsdv';
demand.init(); 

markets = max(demand.marketid);
bootsample = sort(randi(markets, markets, bootreps));
row = (1:size(pk,1))';
bootselection = cell(bootreps, 1);
newmarketid = cell(bootreps, 1);
for r = 1:bootreps
    rows = [];
    newmarket = [];
    for t = 1:markets
        newrows = row(demand.marketid==bootsample(t,r),:);
        rows = [rows; newrows];
        newmarket = [newmarket; ones(size(newrows,1),1) * t];
    end
    bootselection{r} = rows;
    newmarketid{r} = newmarket;
end
pk2 = pk(bootselection{1},:);
bootresults = zeros(length(demand.vars), bootreps);
bootpricech = zeros(6, bootreps); % length(unique(nlMarket.firm))
bootcosts = zeros(6, bootreps); % length(unique(nlMarket.firm))
disp('Market Result')
names = cell(bootreps,1);

parfor_progress(bootreps);
parfor r = 1:bootreps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Nested Logit Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(r)
% Define nested logit demand
newdemand = copy(demand);
newdemand.T = pk(bootselection{r},:);
newdemand.T.newmarketid = newmarketid{r};
newdemand.var.market = 'newmarketid';

newdemand.init(); 
newdemand.options.estimateMethod = 'gmm';
results = newdemand.estimate();
bootresults(:,r) = results{:,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Nested Logit Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


selection = (pk.year==2008 & pk.month==12);
selection = (newdemand.marketid == max(newdemand.marketid ));

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
p0start = nlMarket2.p0;
for j = 1:5
    try
        nlMarket2.equilibrium('Display', 'off');
    catch err
       % disp 'trying again'
        if j == 1
            nlMarket2.p0 = nlMarket.c;
        else
            nlMarket2.p0 = nlMarket.c + (p0start - nlMarket.c)*rand;
        end
    end
end
%market2.fixedPoint(1000)
nlMergerResult = nlMarket.compare(nlMarket2.p);
if r == 1
    names{r} = nlMergerResult.Properties.RowNames;
end
bootpricech(:,r) = nlMergerResult{:,'PriceCh'};
bootcosts(:,r) = nlMergerResult{:,'Costs'};
if display
    disp(results)
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

disp('Costs');
cresults
disp('Price Change');
prresults

priceCh = array2table(bootpricech');
priceCh.Properties.VariableNames = names{1};

costs = array2table(bootcosts');
costs.Properties.VariableNames = names{1};
if saveXLS
    save fn priceCh costs;
xlswrite(fn, presultarray , 'Results', 'B3')
xlswrite(fn, cresultarray , 'Results', 'H3')
xlswrite(fn, names{1} , 'Results', 'A3')
xlswrite(fn, headings , 'Results', 'B2')
xlswrite(fn, names{1} , 'Results', 'G3')
xlswrite(fn, headings , 'Results', 'H2')

xlswrite(fn, {'Price Change'} , 'Results', 'A1')
xlswrite(fn, {'Costs' }, 'Results', 'G1')


xlswrite(fn, names{1}' , 'Price Change', 'A1')
xlswrite(fn, bootpricech' , 'Price Change', 'A2')
xlswrite(fn, names{1}' , 'Costs', 'A1')
xlswrite(fn, bootcosts' , 'Costs', 'A2')
xlswrite(fn, bootresults' , 'Estimate', 'A2')
xlswrite(fn, demand.vars , 'Estimate', 'A1')
end
end
toc

alpha = bootresults(1,:)';
%hist(alpha,40)
