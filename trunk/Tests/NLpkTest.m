classdef NLpkTest < matlab.unittest.TestCase 
    methods
        function result = merge(obj, ces)
            load painkiller9511main2;
            pk.paracetamol = +(pk.substance =='Paracetamol');
            pk.ibuprofen = +(pk.substance =='Ibuprofen');
            pk.asa = +(pk.substance =='ASA');
            pk.constant = ones(size(pk,1),1);
            
            if ces == true
                pk.Xtablets = pk.Xtablets*10e-7;
                demand = NestedLogitDemand(pk);
                demand.settings.ces = true;
                demand.var.marketsize = 'BL_CES';
                demand.var.price = 'Ptablets';
            else
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
            demand.estimate();        
            selection = (pk.year==2008 & pk.month==12);
            
            market = Market(demand);
            market.var.firm = 'firm';
            market.findCosts(selection);
            market.settings.valueShares = false;
            market.settings.weightedAverages = false;
            
            market2 = copy(market);
            market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK';
            market2.p0 = market.p;
            market2.equilibrium(selection);
            %market2.fixedPoint(1000) 
            result = market.compare(market2);
        end
    end
    
    methods (Test)
        function testCes(testCase)   
            result = testCase.merge(true);
            assert(abs(result{1,'Price2'} - 1.743293  )<10e-3)
            assert(abs(result{1,'PriceCh'} - 0.2112071   )<10e-3)
            display '************** Mergersim CES Test passed ******************'
        end
        function testUnit(testCase)   
            result = testCase.merge(false);
            assert(abs(result{1,'Price2'} - 0.533598  )<10e-3)
            assert(abs(result{1,'PriceCh'} - 0.1059824  )<10e-3)
            display '************** Mergersim Unit Test passed ******************'
        end
    end
end