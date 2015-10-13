%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RC bootstrap with cost and demand calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear

conduct = [0, 0.75];
colname = {'NL', 'RC'; 'CES', 'Unit'; '', 'Coll'};
filename = 'cost_calculations.xlsx';

costresults = [];
for rc = 1:2
    for unitdemand = 1:2
        ces = unitdemand == 1
        estimate = true
        
        input.repetitions = 1;
        input.save = true;
        input.optimalIV = true;
        input.newinstruments = true;
        input.withtime = false;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Demand Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if rc == 2
            if ces
                selectRun = [1, 1]; % CES demand
            else
                selectRun = [2, 1]; % Unit demand
            end
            if estimate
                nonlinear = cell(2,1);
                nonlinear(:) = {'paracetamol fizzy branded constant'};
                nind  = [500; 500];
                cesdemand = {ces; ces};
                quaddraws = [10; 12];

                testdata = table2struct(table(nonlinear, cesdemand, quaddraws, nind));
                input.nocons = false;
                input.fn = 'test1';
                input.packDemand = false;
                input.drawmethod = 'quadrature';
                [job,diary] = batchrun(@pkRC, input, testdata, input.repetitions, ...
                    'Parallel', false, 'select', selectRun, ...
                    'randomstream', 99);
                demand = job{1}.demand;
            else
                if ces
                    load demandCES % Saved simulation 3
                else
                    load demandUnit % Saved simulation 7
                end
            end
        else
            pkNL % Estimate NL demand, uses ces parameter
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Standard Merger Simulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(demand.results.estimate)
        
        load painkillers9511main2new;
        pk.paracetamol = +(pk.substance =='Paracetamol');
        pk.ibuprofen = +(pk.substance =='Ibuprofen');
        pk.asa = +(pk.substance =='ASA');
        pk.constant = ones(size(pk,1),1);
        pk.branded = +(pk.brand =='Alvedon')+(pk.brand =='Ipren')+(pk.brand =='Treo');
        pk.fizzy = +(pk.form =='fizzytablet');
        [~,~,pk.date] = unique(pk(:,{'year','month'}));
        pk.time = pk.date + 419;
        pk.marketing = pk.marketing*10^-6;
        pk.lpacksize = log(pk.packsize);
        pk.ldosage = log(pk.dosage);
        [~,~,pk.brandid]=unique(pk.brand);

        pk(isundefined(pk.firm ), :) = [];
        pk.Xtablets2 = pk.Xtablets*10e-7;
        
        pk.firm(pk.firm == 'Ellem') = 'Meda';
        pk.firm(pk.firm == 'Recip') = 'Meda';
        pk.firm(pk.firm == 'Pfizer') = 'McNeil';
        pk.firmsubst = pk.firm;
        pk.firmsubst(pk.brand == 'Ipren') = 'McNeil-Ibu';
        pk.firmsubst(pk.brand == 'Alindrin') = 'Meda-Ibu';
        
        pk([8264 8313 8364],:) = []; % INCORRECT OBSERVATIONS

        if rc == 1
            newdemand = NestedLogitDemand(pk);
        else
            newdemand = MixedLogitDemand(pk);
        end
        newdemand.settings = demand.settings;
        newdemand.var = demand.var;
        newdemand.results = demand.results;
        newdemand.beta = demand.beta;
        newdemand.init();
        newdemand.initSimulation();

        for cc = 1:2
            market = Market(newdemand);
            market.var.firm = 'firm';
            market.settings.conduct = conduct(cc);
            
            market.init();
            market.findCosts();
            if isempty(costresults)
                costresults = market.T(:,{'year','month','firm','product','substance','brand'});
            end
            costresults{:,['c' colname{1, rc} colname{2,unitdemand} colname{3, cc} ]} = market.c;
        end
    end
end

writetable(costresults, filename );    
