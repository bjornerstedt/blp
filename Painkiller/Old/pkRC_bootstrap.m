%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BLP MERGER SIMULATION WITH PAINKILLER DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametric bootstrap merger simulation

drawBeta = true            % Study effect of individual draws in RC
bootreps = 2
fn = 'boot.xls';

display = false;
saveXLS = false;
merge = true;

load painkiller9511main2;
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');
pk.constant = ones(size(pk,1),1);

newdemand = copy(demand);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Standard Merger Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selection = (pk.year==2008 & pk.month==12);

newdemand.initSimulation( selection);

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

    
% Define mixed logit demand
newdemand = copy(demand);
if drawBeta
    newdemand.beta(1) =  bootpars(r, 1);
    newdemand.rc_sigma = bootpars(r, 2:end);
end

selection = (pk.year==2008 & pk.month==12);

newdemand.initSimulation(selection);

nlMarket = Market(newdemand);
nlMarket.var.firm = 'firm';
nlMarket.init();
nlMarket.findCosts();

nlMarket2 = copy(nlMarket);

nlMarket2.firm(nlMarket2.firm == 'AstraZeneca' ) = 'GSK'; 
nlMarket2.p0 = nlMarket.p;
nlMarket2.init(); % Initialize ownership matrix again

porig = nlMarket2.p0; % is p0 changed in Market?
equilibriumAttempts = 5;
for j = 1:equilibriumAttempts
    try
        nlMarket2.equilibrium('Display', 'off');
    catch err
        disp 'trying again'
        if j == 1
            nlMarket2.p0 = nlMarket.c;
        else
            nlMarket2.p0 = nlMarket.c + (porig - nlMarket.c)*rand; % NOTE rand makes parallel nonreplicable
        end
    end
end
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

% disp('Starting beta');
% bootpars
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