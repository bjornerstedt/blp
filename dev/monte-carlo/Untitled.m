testdiff = @(x,y)assert(abs(x - y)<10e-4);
approx = @(x,y)assert(abs((x - y)/y)<1e-3);
demand = {NestedLogitDemand, MixedLogitDemand};

for tc = 1:2
    m = SimMarket(demand{tc});
    m.model.ces = true;
    %         m.estDemand.settings.robust = 'false';
    m.model.endog = false;
    m.model.beta = [-4; 1; 5];
    m.model.markets = 100;
    m.model.randproducts = false;
    m.model.optimalIV = false;
    m.init();
    m.calculateDemand
%     m.findCosts(m.simDemand)
    
%     simResults = m.simulateDemand()
%         m.findCosts(m.simDemand)

%     if tc == 0
%         testdiff(simResults{1, 'sh'} , 0.013209)
%         testdiff(simResults{1, 'p'} , 5.0482)
%     else
%         testdiff(simResults{1, 'sh'} , 0.016533)
%         testdiff(simResults{1, 'p'} , 5.0648)
%     end
    
    result = m.estimate()
%     if tc == 0
%         testdiff(result{'p','Coef'} , -0.9996)
%         testdiff(result{'p','Std_err'} , 0.0049173)
%     else
%         testdiff(result{'p','Coef'} , -1.001)
%         testdiff(result{'p','Std_err'} , 0.0050055)
%     end
end