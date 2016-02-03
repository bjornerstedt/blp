clear

testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);
testtrue = @(x,y)assert(abs((x - y)/y)<1e-3);
testtrue1 = @(x,y)assert(abs((x - y)/y)<1e-2);
testtrue2 = @(x,y)assert(abs((x - y)/y)<2e-2, ...
    'Percentage difference is %f', abs((x - y)/y));

m = SimMarket();
model = m.model;
model.endog = false;
model.markets = 100;
model.products = 5;
model.randproducts = false;
model.optimalIV = false;
model.nonlinear = 'x';
model.epsilon_sigma = 1;
model.rc_sigma = [1];
%% Generate shares using MixedLogitDemandMarket

m = SimMarket(model);
m.demand = MixedLogitDemand2;
% m.demand.settings.nind = 100;
m.init();
% results = m.simulateDemand()
results = m.calculateDemand()
demand = m.estDemand;
demand.beta = model.beta;
demand.data = m.data;
demand.var.quantity = 'sh';
demand.init();
edelta = demand.findDelta(demand.rc_sigma);

md = demand.period{1};
log(md.edelta);

delta = log(md.edelta);
edeltax = log(md.findDelta(md.rc_sigma, md.edelta, 1e-14));

%% Jacobians

F = @(x) md.getShares( delta, x);
ndsdr = jacobian(F, md.rc_sigma);
G = @(x) md.getShares( x, md.rc_sigma);
ndsdd = jacobian(G, delta);
H = @(x,y) md.getShares(x, y);

[dsdd, dsdr] = md.shareJacobians(delta, md.rc_sigma);

dsdr
ndsdr
dsdd
ndsdd

%% Share Hessian calculations

[d2sdd2,d2sdr2,d2sddr,d2sdrd] = md.shareHessians(md.vx, delta, md.rc_sigma);
hyy = hessian(G, delta);
% d2sdd2 

hxx = hessian(F, md.rc_sigma);
% d2sdr2

hxy = crosshessian(@(x,y) H(y,x), md.rc_sigma, delta);
% d2sdrd

hyx = crosshessian(H, delta, md.rc_sigma);
% d2sddr

%% Delta Hessian calculations

% delta

% edelta2 = log(md.findDelta(md.rc_sigma, exp(delta), 1e-14))
% log(md.findDelta(md.rc_sigma, edelta2, 1e-14))
deltajac = md.deltaJacobian( md.rc_sigma, md.edelta)

df = @(x) log(md.findDelta( x, exp(delta), 1e-14));
njac = jacobian(df, md.rc_sigma)

% Numeric and analytic 
nhess = hessian(df, md.rc_sigma);

display 'Numerical Hessian'
nhess

display 'Analytical Hessian'
[jac, hess] = md.deltaDerivatives(delta, md.rc_sigma)

testdiff(hess(1,1), nhess(1,1))
display 'Analytic and numerical Hessians coincide'

m.estDemand.init();
m.estDemand.initEstimation();

rcval = 1.2;
df=@(x)m.estDemand.objective(x);
ng = jacobian(df,rcval)
nh = hessian(df,rcval)
[f, g, h] = m.estDemand.objective(rcval)

if true
    result= zeros(100,2);
    for i = 1:100
        theta = 0.5 + i/100;
        f = m.estDemand.objective(theta);
        result(i,:) = [theta,f];
    end
    plot(result(:,1),result(:,2))
end