clear

m = SimMarket;
m.est.nonlin = 'p';
m.model.endog = true;        % Endog with count instruments or no endog

m.model.randproducts = true; % Let the number of products be random
m.model.beta = [-1; 1; 0];
m.model.rc_sigma = .2;
m.model.x = [2,0];
m.model.x_sigma = [.2,1];
% m.est.findCosts = true; 
m.model.c = 0.9;
m.model.markets = 1000;
m.model.products = 10;

m.run

disp '***************'
selection = m.data.marketid == 1;
m.demand.initSimulation(1);
[~,sh] = m.demand.shares(m.data.p(selection));
display Total Shares
ss = sum(sh)';
Mean = mean(ss);
Min = min(ss);
Max = max(ss);
Std_Dev = std(ss);
disp(table(Mean,Std_Dev,Min,Max))
