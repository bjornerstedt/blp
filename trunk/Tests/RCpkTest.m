% Tests multi market merger simulation for both ces and unit demand
% To run use command: runtests RCpkTest
% Single market tests will be against the next commit 129 when replication 
% of previous multi-market simulations work.  
classdef RCpkTest < matlab.unittest.TestCase 
   methods (TestMethodSetup)
        function setSeed(testCase)
            stream1 = RandStream('mt19937ar','Seed',1);
            RandStream.setGlobalStream(stream1);
        end
   end
   methods
       function result = merge(obj, ces, method, onePeriod)
           load painkiller9511main2;
           pk.paracetamol = +(pk.substance =='Paracetamol');
           pk.ibuprofen = +(pk.substance =='Ibuprofen');
           pk.asa = +(pk.substance =='ASA');
           pk.constant = ones(size(pk,1),1);
           pk.Xtablets2 = pk.Xtablets*10e-7;
           demand = RCDemand2(pk);
           if ces
               demand.settings.ces = true;
               demand.var.marketsize = 'BL_CES';
               demand.var.price = 'Ptablets';
               demand.var.quantity = 'Xtablets2';
           else
               demand.var.marketsize = 'BL_Unit';
               demand.var.quantity = 'Xtablets';
               demand.var.price = 'Ptablets_Real';
           end
           demand.var.nonlinear = 'paracetamol';
           demand.var.market = 'time';
           demand.var.panel = 'product';
           demand.var.exog = ['marketing1 sw sm month2 month3 month4 month5 month6 '...
               'month7 month8 month9 month10 month11 month12'];
           demand.var.instruments = 'num numg numf numfg numhg numfgh';
           demand.settings.paneltype = 'lsdv';
           demand.settings.nind = 300;
           
           demand.config.fptolerance1 = 1e-12; % use lower tolerance for first FP iterations
           demand.config.fptolerance2 = 1e-12; % use maximum tolerance for last iterations
           demand.settings.drawmethod = method;
           demand.init();
           demand.estimate();
           demand.settings.optimalIV = true;
           demand.estimate();
            if onePeriod
                selection = (pk.year==2008 & pk.month == 12);
            else
                selection = (pk.year==2008 );
            end
           
           market = Market(demand);
           market.var.firm = 'firm';
           market.settings.valueShares = false;
           market.settings.weightedAverages = false;
           market.findCosts(selection);
           
           market2 = copy(market);
           market2.firm(market2.firm == 'AstraZeneca' ) = 'GSK';
           market2.p0 = market.p;
           market2.equilibrium(selection);
           result = market.compare(market2);
           result
       end
   end
    
    methods (Test)
        function testHalton2(testCase)
            display '***************** testHalton2 **********************'
            result = testCase.merge(false, 'halton2', true);
            assert(abs(result{1,'Price2'} - 0.54208  )<10e-4)
            assert(abs(result{1,'PriceCh'} - 0.12108 )<10e-3)
            
%             result = testCase.merge(true, 'halton2', true);
%             assert(abs(result{1,'Price2'} - 1.5565  ) <10e-4)
%             assert(abs(result{1,'PriceCh'} - 0.075958  )<10e-4)
        end
        function testHaltonMulti(testCase)
            display '***************** testHaltonMulti **********************'
            result = testCase.merge(false, 'halton', false);
            assert(abs(result{1,'Price2'} - 0.52152  )<10e-4)
            assert(abs(result{1,'PriceCh'} - 0.111   )<10e-3)
            
%             result = testCase.merge(true, 'halton', false);
%             assert(abs(result{1,'Price2'} - 1.4945  ) <10e-4)
%             assert(abs(result{1,'PriceCh'} - 0.058481  )<10e-4)
        end
        function testUnitHalton(testCase)
            display '***************** testUnitHalton **********************'
            result = testCase.merge(false, 'halton', true);
            assert(abs(result{1,'Price2'} - 0.54208  )<10e-4)
            assert(abs(result{1,'PriceCh'} - 0.12108 )<10e-3)
        end
        function testUnitHypercube(testCase)
            display '***************** testUnitHypercube **********************'
            result = testCase.merge(false, 'hypercube', true);
            assert(abs(result{1,'Price2'} - 0.54208  )<10e-4)
            assert(abs(result{1,'PriceCh'} - 0.12108 )<10e-3)
        end
        function testUnitHypercubeMulti(testCase)
            display '***************** testUnitHypercubeMulti **********************'
            result = testCase.merge(false, 'hypercube', false);
            assert(abs(result{1,'Price2'} - 0.52152  )<10e-4)
            assert(abs(result{1,'PriceCh'} - 0.111   )<10e-3)
        end
        function testUnitQuadrature(testCase)
            display '***************** testUnitQuadrature **********************'
            result = testCase.merge(false, 'quadrature', true);
            assert(abs(result{1,'Price2'} - 0.54208  )<10e-4)
            assert(abs(result{1,'PriceCh'} - 0.12108 )<10e-3)
        end
    end
end