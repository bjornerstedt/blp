% pkML2 is basic ML estimations, invoked in parallel
% test==0 is to initialize and return the number of tests to dobatch.m
% Returns a structure of results
% Tests first stage
function results = pkML2(test)
global commondraws;
global repcount;
repcount = 1;
if test == 0
    stream1 = RandStream('mt19937ar','Seed', 99);
    RandStream.setGlobalStream(stream1);

    results = 8;
    return 
end

ces = true
optimalIV = true
newinstruments = true

maxFval = 10^-4;
fval = 10^4;

newIndividualDraws = false;
withtimeTrend = false;

% Either loop over optimalIV to satisfy fval < maxFval 
% or try up to maxTries starting points
if true
    optimalIVmaxRuns = 5;
    maxTries = 1;
else
    maxTries = 5;
    optimalIVmaxRuns = 1;
end
if ces
    price = 'Ptablets ';
else
    price = 'Ptablets_Real ';
end
commonRC = 'constant ' % RC in all treatments
demandcase = {['paracetamol ']};
demandcase{2} = ['paracetamol fizzy '];
demandcase{3} = ['paracetamol fizzy branded '];
demandcase{4} = ['paracetamol fizzy ', price ];
demandcase{5} = ['paracetamol ',commonRC];
demandcase{6} = ['paracetamol fizzy ',commonRC];
demandcase{7} = ['paracetamol fizzy branded ',commonRC];
demandcase{8} = ['paracetamol fizzy ', price ,commonRC];
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

% ***** Loop for completely new starting points *****
tryit=0;
while fval > maxFval & tryit <= maxTries

if ces
    pk.Xtablets = pk.Xtablets*10e-7;
    demand = CesMixedLogitDemand(pk);
    demand.marketsize = 'BL_CES';
    demand.price = 'Ptablets'; 
else
    demand = MixedLogitDemand(pk);
    demand.marketsize = 'BL_Unit';
    demand.price = 'Ptablets_Real'; 
end

demand.nonlinear = demandcase{test};

demand.quantity = 'Xtablets';
demand.market = 'time';
demand.panel = 'product';

%demand.panel = 'brandid';
demand.exog = ['marketing sw sm month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
if withtimeTrend
    demand.exog = ['time ' demand.exog];
end
if newinstruments
    demand.instruments = ['i1_con i2_con i1_ldosage i2_ldosage i1_lpacksize '...
        'i2_lpacksize i1_form2 i2_form2 i1_substance2 i2_substance2 '...
        'i1_substance3 i2_substance3'];
else
    demand.instruments = 'num numg numf numfg numhg numfgh';
end

demand.settings.drawmethod = 'hypercube';
demand.settings.marketdraws = true;
demand.settings.nind = 500;
demand.paneltype = 'lsdv';
demand.settings.quaddraws = 5;

% ***** Loop for up to optimalIVmaxRuns repeated optimalIV iterations
ivit=0;
ivflag = optimalIV; 
firsttries = 1;
while fval > maxFval & ivit <= optimalIVmaxRuns | ivflag
    if ivit==0
        disp('****** First round estimation ******');
        % Save draws of individuals in test 1 to global
        if test == 1 | newIndividualDraws | repcount == 1
            demand.rc_sigma = [];
            demand.init();
            commondraws = demand.draws;
            repcount = repcount + 1;
        else
            demand.draws = commondraws;
            demand.init();
        end
        results.estimate = demand.estimate();
        results.estimate2 = demand.estimate(); % Hack: If no second stage
        results.fval = demand.results.fval;
        if demand.results.cond > 10^8 & firsttries <= 5
            disp(['****** Repeating estimations as cond=',num2str(demand.results.cond), ' *******']);
            firsttries = firsttries + 1;
        else
            ivit = ivit+1;
        end
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
            if demand.results.isreal
                fval = demand.results.fval;
            end
        end
        ivit = ivit+1;
    end
end
tryit = tryit + 1;
end
results.fval = [results.fval; fval; ivit-1; demand.results.cond; ...
    demand.results.isreal; demand.results.sdreg];
merger = Merger();
merger.selection = pk.year==2008 & pk.month==12;
merger.firmvar = 'firm';
merger.buyer = 'GSK';
merger.seller = 'AstraZeneca';
results.merger = merger.merge(demand);

delete(demand);
%delete(merger);

end


% while fval > maxFval & i < 10
%     if i==0
%         % Save draws of individuals in test 1 to global
%         if test == 1 | newIndividualDraws 
%             demand.init();
%             commondraws = demand.draws;
%         else
%             demand.draws = commondraws; 
%             demand.init();
%         end
%     else
%         disp(['**** Repeating estimations as fval=',num2str(fval)]);
%     end
%     results.estimate = demand.estimate();
%     results.fval = demand.results.fval;
% %    results.merger = merger.merge(demand);
%     if optimalIV & i == 0
%         demand.settings.optimalIV = true;
%         results.estimate2 = demand.estimate();
%     end
%     results.merger = merger.merge(demand);
%     i = i+1;
%     fval = demand.results.fval;
%     results.fval = [results.fval; fval];
% end
