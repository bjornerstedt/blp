%% Calculating the delta Hessian array
% 
% This calculation tests the derivation of the delta Hessian array from 
% share Jacobians and Hessians as described in the text:
%    multi-dimensional derivatives.lyx

%% Definitions

% clear all
testdiff = @(x,y)assert(abs(x - y)<10e-4, 'Obtained %f instead of %f' ,x, y);

import HessianCalc.*

stream1 = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(stream1);
rc = 2;
pr = 4;
n = 10;
v = randn(1, n , rc);
z = abs(randn(pr, 1));
vx = bsxfun(@times, z, v);

delta = randn(pr,1);
rc_sigma = randn(rc,1);
[s, S] = share(vx, delta, rc_sigma);

%% Jacobians

F = @(x) share(vx, delta, x);
jac1 = jacobian(F, rc_sigma);

G = @(x) share(vx, x, rc_sigma);
H = @(x,y) share(vx, x, y);
jac2 = jacobian(G, delta);
jac1
jac2

[dsdd, dsdr] = shareJacobians(vx, delta, rc_sigma)
% sum(dS_dRhoCalc(S, vx), 2)
% Note that dividing by obj.settings.nind cancels out.
dddr = -dsdd\dsdr;

%% Share Hessian calculations

[d2sdd2,d2sdr2,d2sddr,d2sdrd] = share_hessians(S, vx);
hyy = hessian(G, delta)
d2sdd2 

hxx = hessian(F, rc_sigma)
d2sdr2

hxy = crosshessian(@(x,y) H(y,x), rc_sigma, delta)
d2sdrd

hyx = crosshessian(H, delta, rc_sigma)
d2sddr

%% Delta Hessian calculations

delta
log(finddelta(s, vx, rc_sigma, delta))

df = @(x) log(finddelta(s, vx, x, delta));

% Numeric and analytic jacobians give the same result

ddelta_dsigma = jacobian(df, rc_sigma);

% Numeric and analytic 
nhess = hessian(df, rc_sigma);

id = arraymult(arraymult(d2sdd2, dddr, 2), dddr);
cd1 = arraymult(d2sdrd, dddr,2);
cd2 = arraymult(d2sddr, dddr);

invdsdd = inv(dsdd);
analder = -(id + cd1 + cd2 + d2sdr2);

% The analytic and numerical Hessians are almost identical!
analhess = arraymult(invdsdd, analder)

display 'Numerical Hessian'
nhess

display 'Analytical Hessian'
[~, hess] = deltaDerivatives(vx, delta, rc_sigma)

testdiff(hess(1,1,1), 0.017144)
testdiff(hess(1,2,2), 0.2657)
display 'Tests passed'