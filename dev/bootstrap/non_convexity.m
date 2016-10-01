%% Merger with estimated demand
% display '**********************   Estimated  *************************'
% m = SimMarket();
% m.model.markets = 50;
% m.demand = NLDemand;
% m.demand.alpha = 2;
% m.model.endog = true;
% m.model.randomProducts = true;
% m.model.firm = [1,1,2,2,3];
% 
% % pl = zeros(10,2);
% % for i = 1:10;
% %     j =  1 * i;
% %     pl(i,:) = [j, m.calibrate(.5, j)];
% % end
% % plot(pl(:,1),pl(:,2))
% 
% [mk, beta0] = m.calibrate(.5)
% sum(mk{:,2})
% 
% m.demand.data.firm2 = m.data.firm;
% m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;
% 
% display(m.model)
% results = m.demand.estimate()
% 
% 
% m.findCosts();
% market = m.market;
% 
% market.summary()
% market.summary('Selection', m.data.marketid == 1)
% 
% 
% market2 = copy(market);
% market2.var.firm = 'firm2';
% market2.equilibrium();
% 
% mergerResult1 = summary(market, market2)

%% Calculate the true price increase
display '**********************   TRUE  *************************'

res = zeros(20, 2);
for i = 1:40
    alpha = 1/20 * i;
    m = SimMarket();
    m.model.markets = 2;
    m.demand = NLDemand;
    m.demand.alpha = alpha;
%     m.model.endog = true;
%     m.model.randomProducts = true;
    m.model.firm = [1,1,2,2,3];
    m.model.eta = 0;
    m.model.xi = 0;
    m.model.x_vcv = [0 0];
    
    m.calibrate(.5);
    
    m.demand.data.firm2 = m.data.firm;
    m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;
    
    % display(m.model)
    
    m.findCosts();
    market = m.market;
    
    market2 = copy(market);
    market2.var.firm = 'firm2';
    market2.equilibrium();
    
    mergerResult = summary(market, market2);
    res(i, 1) = alpha;
    res(i, 2) = mergerResult{1, 'PriceCh'};
end
plot(res(:,1),res(:,2))

%% Aggregate shocks

bdev = zeros(20, 2);
for i = 1:40
    alpha = .5;
    m = SimMarket();
    m.model.markets = 5;
    m.demand = NLDemand;
    m.demand.alpha = alpha;
    m.model.firm = [1,1,2,2,3];
    m.model.eta = 0;
    m.model.xi = 0;
    m.model.x_vcv = [0 0];
    
    m.varyIntercept(i/4);
    
    m.demand.data.firm2 = m.data.firm;
    m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;
    
    % display(m.model)
    
    m.findCosts();
    market = m.market;
    
    market2 = copy(market);
    market2.var.firm = 'firm2';
    market2.equilibrium();
    
    mergerResult = summary(market, market2);
    bdev(i, 1) = i/4;
    bdev(i, 2) = mergerResult{1, 'PriceCh'};
end
plot(bdev(:,1),bdev(:,2))