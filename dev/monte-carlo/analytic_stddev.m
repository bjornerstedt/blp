m = SimMarket(NestedLogitDemand);
m.model.markets = 1000;
m.model.randproducts = false;
m.init
X = data{:,{'p' , 'x' , 'constant'}};
mean(X)
1000*inv(X'*X)