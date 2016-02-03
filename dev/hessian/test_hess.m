%% Testing the second order implicit function theorem
% 
% The purpose of this script is to test the derivation of the 
% three-dimensional Hessian array. Here simple functions are used to check
% the validity of the derivatives.
% 
%% Real valued function
%
% We begin by analysis of the real valued function $f(x,y)$ 
% 
% $$f(x, y) = x^2*y^2 - 1 = 0$$
% 
% The implicitly defined function is $y = 1/x$ with first derivative given
% by $-x^{-2}$ and second derivative $2x^{-3}$

t = 3;
x0 = t;
y0 = 1/t;
f = @(x, y) x^2*y^2;
d2fdx2 = 2*y0^2;
d2fdy2 = 2*x0^2;
d2fdxy = 4*x0*y0;

truejac = -x0^-2;
truehess = 2*x0^-3;

f(x0, y0);
dfx = jacobian( @(x) f(x, y0), x0);
dfy = jacobian( @(y) f( x0, y), y0);

truejac
dydx = -dfy\dfx

hx1 = jacobian( @(w)jacobian( @(x) f(x, y0), w), x0);
hxy1 = jacobian( @(y)jacobian( @(x) f(x, y), x0), y0);
hyx1 = jacobian( @(x)jacobian( @(y) f(x, y), y0), x0);
hy1 = jacobian( @(w)jacobian( @(y) f(x0, y), w), y0);

hh1 = hy1*dydx*dydx;
hh2 = hyx1*dydx;
% Hessian
ders1 = hy1*dydx*dydx + hyx1*dydx + hxy1*dydx + hx1;

truehess
d2ydx2 = -ders1/dfy

%% Simple vector valued function
f = @(x, y) [x(1)^2*y(1)^2; x(2)^2*y(2)^2];

x0 = [t;t];
y0 = [1/t;1/t];

f(x0, y0);
dfx = jacobian( @(x) f(x, y0), x0);
dfy = jacobian( @(y) f( x0, y), y0);

truejac
dydx = -dfy\dfx

hx = hessian( @(x) f(x, y0), x0); 
hy = hessian(@(y) f( x0, y), y0);
hxy = crosshessian(f, x0, y0);
hyx = crosshessian(@(x,y) f(y,x), y0, x0);

ders = arraymult(arraymult(hy, dydx, 2), dydx) + ...
arraymult(hxy, dydx) + ...
arraymult(hyx, dydx, 2) + ...
hx;

% The hessian is the same

truehess
ah = arraymult( -inv(dfy), ders)

%% Vector valued function with interactions
%
% Now we model 
%
% $$y_1 = \frac{1}{2x_1 + x_2} $$
%
% $$y_2 = \frac{1}{x_1 + 2x_2} $$
%
% using the function
%
% $$f(x,y) = [(2x_1 + x_2)^2y_1^2, (x_1 + 2x_2)^2y_2^2]-[1,1] = [0,0]$$
%

x0 = [1;1];
y0 = [1/3;1/3];

% The implicitly defined function has derivatives given by 

truejac = -[2,1;1,2]*(2*x0(1)+x0(2))^-2
truehess = 2*[4,2;1,2]*(2*x0(1)+x0(2))^-3

f = @(x, y) [(2*x(1)+x(2))^2*y(1)^2; (x(1)+2*x(2))^2*y(2)^2];

f(x0, y0)

dfx = jacobian( @(x) f(x, y0), x0)
dfy = jacobian( @(y) f( x0, y), y0)

truejac
dydx = -dfy\dfx

hx = hessian( @(x) f(x, y0), x0) 
hy = hessian(@(y) f( x0, y), y0)
hxy = crosshessian(f, x0, y0)
hyx = crosshessian(@(x,y) f(y,x), y0, x0)

hh1 + hh2
arraymult(arraymult(hy, dydx, 2), dydx) +arraymult(hyx, dydx) 

ders1
ders = arraymult(arraymult(hy, dydx, 2), dydx) + ...
arraymult(hyx, dydx) + ...
arraymult(hxy, dydx, 2) + ...
hx

% The hessian is the same

truehess
ah = arraymult( -inv(dfy), ders)

%% Vector valued function with x of lower dimensionality
%
% Now we model 2x3 -> 3 dimensions
%
% $$y_1 = \frac{1}{2x_1 + x_2} $$
%
% $$y_2 = \frac{1}{x_1 + 2x_2} $$
%
% $$y_3 = \frac{1}{x_1} $$
%
% using the function
%
% $$f(x,y) = [(2x_1 + x_2)^2y_1^2, (x_1 + 2x_2)^2y_2^2, x_1^2 y_3^2]-[1,1,1] = [0,0,0]$$
%
t=2;
x0 = [t;t];
y0 = [1/(3*t);1/(3*t);1/t];

f = @(x, y) [(2*x(1)+x(2))^2*y(1)^2; (x(1)+2*x(2))^2*y(2)^2; x(1)^2*y(3)^2 ];

f(x0, y0)

% The implicitly defined function has derivatives given by 

truejac = -[[2,1;1,2]*(2*x0(1)+x0(2))^-2; [1,0]*x0(1)^-2]
truehess = 2*[4,2;1,2]*(2*x0(1)+x0(2))^-3


dfx = jacobian( @(x) f(x, y0), x0)
dfy = jacobian( @(y) f( x0, y), y0)

truejac
dydx = -dfy\dfx

hx = hessian( @(x) f(x, y0), x0) 
hy = hessian(@(y) f( x0, y), y0)
hxy = crosshessian(f, x0, y0)

hyx = crosshessian(@(x,y) f(y,x), y0, x0)

arraymult(arraymult(hy, dydx, 2), dydx) 

hh1 + hh2
arraymult(arraymult(hy, dydx, 2), dydx) + arraymult(hyx, dydx) 

ders1
ders = arraymult(arraymult(hy, dydx, 2), dydx) + ...
arraymult(hyx, dydx) + ...
arraymult(hxy, dydx, 2) + ...
hx

% The hessian is the same

truehess
ah = arraymult( -inv(dfy), ders)
