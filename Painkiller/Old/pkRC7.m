function results = pkRC(input, testdata, repetition)
% pkRC implements a parallel estimation
% repetition==0 is to initialize and return a struct to include in
% testdata.
% Returns a structure of results

optimalIV = true
newinstruments = true

maxFval = 10^-4;
fval = 10^4;

newIndividualDraws = false;
withtime = false;

% disp(['*** Test ', num2str(test), ', Non-linear: ', demandcase{test}, ' ***'])

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
optimalIVmaxRuns = 5;

if testdata.ces
    pk.Xtablets = pk.Xtablets*10e-7;
    demand = CesMixedLogitDemand(pk);
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
else
    demand = MixedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
end

demand.var.nonlinear = testdata.nonlinear;

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
demand.settings.marketdraws = input.marketdraws;
demand.settings.nind = 1000;
demand.settings.paneltype = 'lsdv';
demand.settings.quaddraws = 5;

if repetition == 0
 %   disp(['************ Executing test  ',num2str(test), ' *************']);
    stream1 = RandStream('mt19937ar','Seed', 99);
    RandStream.setGlobalStream(stream1);
    demand.init();
    results = testdata;
    results.commondraws = demand.draws;
    results.rc_sigma = randn(input.repetitions, size(demand.rc_sigma, 1));
    return
else
    if testdata.commonDraws
        demand.draws = testdata.commondraws;
    end
    demand.rc_sigma = testdata.rc_sigma(repetition,:);
    demand.init();
end
disp(['****** Repetition ',num2str(repetition), ' *******']);
results.estimate = demand.estimate();
results.fval = demand.results.fval;
% ***** Loop for up to optimalIVmaxRuns repeated optimalIV iterations
ivit=1;
while fval > maxFval & ivit <= optimalIVmaxRuns 
    if ivit == 1
        disp('****** Optimal IV estimation ******');
    else
        disp(['****** Repeating estimations as fval=',num2str(fval), ' *******']);
    end
    demand.settings.optimalIV = true;
    if ivit > 1
        demand.rc_sigma = randn(size(demand.rc_sigma));
    end
    results.estimate2 = demand.estimate();
    fval = demand.results.fval;
    ivit = ivit+1;
end
results.fval = [results.fval; fval; ivit-1; demand.results.cond; ...
    demand.results.sdreg];
merger = Merger();
merger.selection = pk.year==2008 & pk.month==12;
merger.var.firm = 'firm';
merger.buyer = 'GSK';
merger.seller = 'AstraZeneca';
results.merger = merger.merge(demand);

delete(demand);
%delete(merger);

end
