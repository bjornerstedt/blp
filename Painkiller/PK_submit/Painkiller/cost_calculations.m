%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Painkiller NL and RC cost calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);

filename = 'cost_calculations2.csv';
estimate = true;

load painkillers;

costresults = pk(:,{'year','month','firm','product','substance','brand'});
for nl_rc = 1:2
    if nl_rc == 2
        disp '**************************************************************'
        disp '******************  Mixed Logit Demand  **********************'
        disp '**************************************************************'
        disp ' '
        if estimate
            demand = MixedLogitDemand(pk(pk.year<=2008, :));
            display 'CES Demand';
            demand.settings.ces = true;
            demand.var.marketsize = 'BL_CES';
            demand.var.price = 'Ptablets';
            demand.var.quantity = 'Xtablets1';
            demand.var.market = 'time';
            demand.var.panel = 'product';
            
            demand.var.nonlinear = 'paracetamol fizzy branded constant';
%             demand.var.nonlinear = 'paracetamol';
            demand.var.exog = ['marketing sw sm month2 month3 month4 month5 month6 '...
                'month7 month8 month9 month10 month11 month12'];
            demand.var.instruments = ['i1_con i2_con i1_ldosage i2_ldosage '...
                'i1_lpacksize i2_lpacksize i1_form2 i2_form2 i1_substance2 '...
                'i2_substance2 i1_substance3 i2_substance3'];
            demand.settings.drawmethod = 'quadrature';
            demand.settings.paneltype = 'lsdv';
            demand.settings.quaddraws = 10;
            
            demand.init();
            demand.estimate();
            demand.settings.optimalIV = true;
            demand.estimate();
            save /Users/jonasbjornerstedt/Documents/demandCES.mat demand;
        else % Saved estimation
            load /Users/jonasbjornerstedt/Documents/demandCES.mat
            disp(demand.results.estimate)
        end
    else
        disp '**************************************************************'
        disp '*****************  Nested Logit Demand  **********************'
        disp '**************************************************************'
        disp ' '
        demand = NestedLogitDemand(pk(pk.year<2009, :));
        demand.settings.ces = true;
        demand.var.nests = 'form substance';
        demand.var.exog = ['marketing sw sm time month2 month3 month4 month5 month6 '...
            'month7 month8 month9 month10 month11 month12'];
        disp 'CES Demand'
        demand.var.marketsize = 'BL_CES';
        demand.var.price = 'Ptablets';
        demand.var.quantity = 'Xtablets1';
        demand.var.market = 'date';
        demand.var.panel = 'product';
        demand.var.instruments = 'num numg numf numfg numhg numfgh';
        demand.settings.paneltype = 'lsdv';
        demand.init();
        demand.estimate()
    end
    demand.T = pk; % Replace data with whole time period
    demand.init();
    demand.initSimulation();
    conduct = [0, 0.75];
    colname = {'NL', 'RC';    '', 'Coll'};
    for cc = 1:2
        market = Market(demand);
        market.var.firm = 'firm';
        market.settings.conduct = conduct(cc);
        market.init();
        market.findCosts();
        costresults{:,['c' colname{1, nl_rc} 'CES' colname{2, cc} ]} = market.c;
    end
end
% testval = mean(costresults{costresults.year == 2008 & costresults.month == 12,8:9});
% assert(abs(testval(1)-0.386430637760993)<1e-7)
% assert(abs(testval(2)-0.880465843450507)<1e-7)
% display '**************** TESTS PASSED *********************'
writetable(costresults, filename );
