%% Merger with estimated demand
if false
display '**********************   Estimated  *************************'
m = SimMarket();
m.model.markets = 50;
m.demand = NLDemand;
m.demand.alpha = .5;
m.model.endog = true;
m.model.randomProducts = true;
m.model.firm = [1,1,2,2,3];

% pl = zeros(10,2);
% for i = 1:10;
%     j =  1 * i;
%     pl(i,:) = [j, m.calibrate(.5, j)];
% end
% plot(pl(:,1),pl(:,2))

[mk, beta0] = m.calibrate(.5)
sum(mk{:,2})

m.demand.data.firm2 = m.data.firm;
m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;

display(m.model)
results = m.demand.estimate()


m.findCosts();
market = m.market;

market.summary()
market.summary('Selection', m.data.marketid == 1)


market2 = copy(market);
market2.var.firm = 'firm2';
market2.equilibrium();

mergerResult1 = summary(market, market2)
end
%% Calculate the true price increase
display '**********************   Vary alpha and beta0  *************************'

res = zeros(20, 2);
alpha = linspace(.3, 1, 40);
alpha = .5 + .15 * randn( 40, 1);
beta0Vec = linspace(-1, 5, 40);
qVec = linspace(.1, .9, 40);
for i = 1:40
    m = SimMarket();
    m.model.markets = 1;
    m.demand = NLDemand;
    m.demand.alpha = .5;
    m.demand.alpha = alpha(i);
    m.model.firm = [1,1,2,2,3];
    m.model.eta = 0;
    m.model.xi = 0;
    m.model.x_vcv = [0 0];

    if false
          mk = m.calibrate(.5);
%         mk = m.calibrate(qVec(i));
    else
        m.model.beta(2) = beta0Vec(i);
        m.model.beta(2) = 1;
        mk = m.create();
    end
    
    m.demand.data.firm2 = m.data.firm;
    m.demand.data.firm2(m.demand.data.firm2 == 2 ) = 1;
    
    % display(m.model)
    
    m.findCosts();
    market = m.market;
    
    market2 = copy(market);
    market2.var.firm = 'firm2';
    market2.equilibrium();
    
    mergerResult = summary(market, market2);
    res(i, 1) = alpha(i);
    res(i, 2) = m.model.beta(2);
    res(i, 3) = mergerResult{1, 'PriceCh'};
    res(i, 4) = sum(mk{:, 'q'});
    res(i, 5) = qVec(i);
end
mean(res(:,3))
% plot(res(:,1),res(:,3))
% plot(res(:,4),res(:,3))
return
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