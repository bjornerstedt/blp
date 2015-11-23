function results = pkRC(input, testdata, repetition, randstream)
% pkRC implements a parallel estimation
% repetition==0 is to initialize and return a struct to include in
% testdata.
% Returns a structure of results

optimalIV =  true;
newinstruments = testdata.newinstruments;

maxFval = 10^-4;
fval = 10^4;

newIndividualDraws = false;
withtime = testdata.withtime;

% disp(['*** Test ', num2str(test), ', Non-linear: ', demandcase{test}, ' ***'])

load painkillers9511main2new;
%load painkiller9511main2;
pk.Xtablets2 = pk.Xtablets*10e-7;

pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
pk.fizzy = +(pk.form =='fizzytablet');
[~,~,pk.date] = unique(pk(:,{'year','month'}));
%pk.time = pk.date + 419;
pk(pk.year>2008, :) = [];
pk.marketing = pk.marketing*10^-6;
pk.lpacksize = log(pk.packsize);
pk.ldosage = log(pk.dosage);
[~,~,pk.brandid]=unique(pk.brand);

pk.firm(pk.firm == 'Ellem') = 'Meda';
pk.firm(pk.firm == 'Recip') = 'Meda';
pk.firm(pk.firm == 'Pfizer') = 'McNeil';
pk.firmsubst = pk.firm;
pk.firmsubst(pk.brand == 'Ipren') = 'McNeil-Ibu';
pk.firmsubst(pk.brand == 'Alindrin') = 'Meda-Ibu';

% ************** Loop for completely new starting points ******************
% Either loop over optimalIV to satisfy fval < maxFval 
% or try up to maxTries starting points
optimalIVmaxRuns = 1;

if testdata.cesdemand
    demand = MixedLogitDemand(pk);
    demand.settings.ces = true;
    demand.var.marketsize = 'BL_CES';
    demand.var.price = 'Ptablets'; 
    demand.var.quantity = 'Xtablets2';
else
    demand = MixedLogitDemand(pk);
    demand.var.marketsize = 'BL_Unit';
    demand.var.price = 'Ptablets_Real'; 
    demand.var.quantity = 'Xtablets';
end

demand.var.nonlinear = testdata.nonlinear;

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
demand.settings.paneltype = 'lsdv';
demand.settings.marketdraws = true;
demand.settings.nind = testdata.nind;
demand.settings.quaddraws = 5;
demand.config.randstream = randstream;
demand.settings.nocons = input.nocons;

if repetition == 0
 %   disp(['************ Executing test  ',num2str(test), ' *************']);
%     stream1 = RandStream('mt19937ar','Seed', 99);
%     RandStream.setGlobalStream(stream1);
    demand.init();
    results = testdata; % returned results replaces testdata
%     results.draws = demand.draws; % Save draws in testdata
    return
else
%     if testdata.commonDraws
%         demand.draws = testdata.draws;
%     end
%    demand.rc_sigma = testdata.rc_sigma(repetition,:);
    demand.init();
end

disp(['****** Repetition ',num2str(repetition), ' *******']);
demand.estimate();
results.estimate(1) =  demand.results;

if input.packDemand
    results.demand(1) = demand.pack();
else
    results.demand = copy(demand);
end

% ***** Loop for up to optimalIVmaxRuns repeated optimalIV iterations
if optimalIV
    ivit=1;
    while fval > maxFval & ivit <= optimalIVmaxRuns
        if ivit == 1
            disp('****** Optimal IV estimation ******');
        else
            disp(['****** Repeating estimations as fval=',num2str(fval), ' *******']);
        end
        demand.settings.optimalIV = true;
        if ivit > 1
            demand.rc_sigma = randstream.randn(size(demand.rc_sigma));
        end
        demand.estimate();
        results.estimate(2) =  demand.results;
        fval = demand.results.other.fval;
        ivit = ivit+1;
    end
    if input.packDemand
        results.demand(2) = demand.pack();
    else
        results.demand = copy(demand);
    end
end

merger = Merger();
merger.selection = pk.year==2008 & pk.month==12;
merger.var.firm = 'firm';
merger.buyer = 'GSK';
merger.seller = 'AstraZeneca';
results.settings = demand.settings;
results.var = demand.var;
results.merger = merger.merge(demand);
results.equilibrium = merger.results.isequilibrium;
if ~input.packDemand
    delete(demand);
end
delete(merger);

end
