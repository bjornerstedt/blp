% painkillerMLpref is used to invoke tests in parallel
% test==0 is to initialize and return the number of tests to dobatch.m
% Returns a structure of results
function results = pkML(test)
global commondraws;
global repcount;
repcount = 1;
if test == 0
    stream1 = RandStream('mt19937ar','Seed', 99);
    RandStream.setGlobalStream(stream1);

    results = 11;
    return 
end

ces = true
optimalIV = true
newinstruments = true

maxFval = 10^-4;
fval = 10^4;

newIndividualDraws = false;
withtime = false;

if ces
    price = 'Ptablets ';
else
    price = 'Ptablets_Real ';
end
commonRC = 'constant ' % RC in all treatments
demandcase = {['paracetamol ibuprofen ',commonRC]};
demandcase{2} = ['paracetamol ibuprofen asa ',commonRC];
demandcase{3} = ['paracetamol ibuprofen asa fizzy ',commonRC];
demandcase{4} = ['paracetamol ibuprofen fizzy branded ',commonRC];
demandcase{5} = ['paracetamol ibuprofen fizzy lpacksize ',commonRC];
demandcase{6} = ['paracetamol ibuprofen fizzy ldosage ',commonRC];
demandcase{7} = ['paracetamol ibuprofen fizzy ', price ,commonRC];
demandcase{8} = ['paracetamol fizzy branded lpacksize ',commonRC];
demandcase{9} = ['paracetamol fizzy lpacksize ldosage ',commonRC];
demandcase{10} = ['paracetamol fizzy branded ldosage ',commonRC];
demandcase{11} = ['paracetamol fizzy ldosage ', price ,commonRC];
disp(['*** Test ', num2str(test), ', Non-linear: ', demandcase{test}, ' ***'])

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


% ************** Loop for completely new starting points ******************
% Either loop over optimalIV to satisfy fval < maxFval 
% or try up to maxTries starting points
if true
    optimalIVmaxRuns = 5;
    maxTries = 1;
else
    maxTries = 5;
    optimalIVmaxRuns = 1;
end

tryit=0;
while fval > maxFval & tryit <= maxTries

if ces
    pk.Xtablets = pk.Xtablets*10e-7;
    demand = CesMixedLogitDemand(pk);
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
else
    demand = MixedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
end

demand.var.nonlinear = demandcase{test};

demand.var.quantity = 'Xtablets';
demand.var.market = 'time';
demand.var.panel = 'product';

%demand.var.panel = 'brandid';
demand.var.exog = ['marketing sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
if withtime
    demand.var.exog = ['time ' demand.var.exog];
end
if newinstruments
    demand.var.instruments = ['i1_con i2_con i1_ldosage i2_ldosage i1_lpacksize '...
        'i2_lpacksize i1_form2 i2_form2 i1_substance2 i2_substance2 '...
        'i1_substance3 i2_substance3'];
else
    demand.var.instruments = 'num numg numf numfg numhg numfgh';
end

demand.settings.drawmethod = 'hypercube';
demand.settings.marketdraws = false;
demand.settings.nind = 500;
demand.settings.paneltype = 'lsdv';
demand.settings.quaddraws = 5;

% ***** Loop for up to optimalIVmaxRuns repeated optimalIV iterations
ivit=0;
ivflag = optimalIV; 
while fval > maxFval & ivit <= optimalIVmaxRuns | ivflag
    if ivit==0
        % Save draws of individuals in test 1 to global
        if test == 1 | newIndividualDraws | repcount == 1
%            demand.rc_sigma = [0.094622 3.7391 0.89851]';
            demand.init();
            commondraws = demand.draws;
            repcount = repcount + 1;
        else
            demand.draws = commondraws;
            demand.init();
        end
        results.estimate = demand.estimate();
        results.fval = demand.results.fval;
    else
        if optimalIV
            if ivit == 1
                disp('****** Optimal IV estimation ******');
            else
                disp(['****** Repeating estimations as fval=',num2str(fval), ' *******']);
            end
            demand.settings.optimalIV = true;
            ivflag = false; % Continue iterations only if fval > maxFval
            results.estimate2 = demand.estimate();
            fval = demand.results.fval;
        end
    end
    ivit = ivit+1;
end
tryit = tryit + 1;
end
results.fval = [results.fval; fval; ivit];
merger = Merger();
merger.selection = pk.year==2008 & pk.month==12;
merger.var.firm = 'firm';
merger.buyer = 'GSK';
merger.seller = 'AstraZeneca';
results.merger = merger.merge(demand);

delete(demand);
%delete(merger);

end
