% Quick RC test
ces = false
method = 'hypercube'
testOnePeriod = true
optimalIV = false

           load painkiller9511main2;
           pk.paracetamol = +(pk.substance =='Paracetamol');
           pk.ibuprofen = +(pk.substance =='Ibuprofen');
           pk.asa = +(pk.substance =='ASA');
           pk.constant = ones(size(pk,1),1);
           if ces
               pk.Xtablets = pk.Xtablets*10e-7;
               demand = MixedLogitDemand(pk);
               demand.settings.ces = true;
               demand.var.marketsize = 'BL_CES';
               demand.var.price = 'Ptablets';
           else
               demand = MixedLogitDemand(pk);
               demand.var.marketsize = 'BL_Unit';
               demand.var.price = 'Ptablets_Real';
           end
           demand.var.nonlinear = 'paracetamol';
           demand.var.quantity = 'Xtablets';
           demand.var.market = 'time';
           demand.var.panel = 'product';
           demand.var.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
               'month7 month8 month9 month10 month11 month12'];
           demand.var.instruments = 'num numg numf numfg numhg numfgh';
           demand.settings.paneltype = 'lsdv';
           demand.settings.nind = 300;
           demand.config.test = true; % Run hypercube/Halton with weights instead of 1/N
           
           demand.settings.fptolerance1 = 1e-12; % use lower tolerance for first FP iterations
           demand.settings.fptolerance2 = 1e-12; % use maximum tolerance for last iterations
           demand.settings.drawmethod = method;
demand.rc_sigma = 4.66;
           demand.init();
           demand.estimate();
           if optimalIV
           demand.settings.optimalIV = true;
           demand.estimate();
           end
            if testOnePeriod
                selection = (pk.year==2008 & pk.month == 12);
            else
                selection = (pk.year==2008 );
            end
           
           market = Market(demand);
           market.var.firm = 'firm';
           market.findCosts(selection);
           
           market2 = copy(market);
           market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK';
           market2.p0 = market.p;
           market2.equilibrium(selection);
           result = market.compare(market2.p);
           result
           pvect = [pk.year, pk.month, market2.p];
           pvect = pvect(selection,:);

if testOnePeriod
    if optimalIV
        assert(abs(result{1,'Price2'} - 0.54208  )<10e-4)
        assert(abs(result{1,'PriceCh'} - 0.12108 )<10e-3)
    else
        assert(abs(result{1,'Price2'} - 0.67219   )<10e-4)
        assert(abs(result{1,'PriceCh'} - 0.39159 )<10e-3)
    end
else
    assert(abs(result{1,'Price2'} - 0.52152  )<10e-4)
    assert(abs(result{1,'PriceCh'} - 0.111   )<10e-3)
end
            display '*********** Mergersim Unit Demand Test passed ****************'
            
