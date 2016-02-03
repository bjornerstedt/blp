%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding Nested Logit FOC Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% The calculations of the FOC Jacobian is done on the last period of the 
% Painkillers dataset.
% Put shareHessian in NestedLogitDemand, 
%   Make compatible with MixedLogitDemand.shareHessian
% Put foc_jacobian in Market

%%  Estimate demand and calculate market costs
clear;

load 'painkillers9511main2new.mat';
getpainkiller
pk.paracetamol = +(pk.substance =='Paracetamol');
pk.ibuprofen = +(pk.substance =='Ibuprofen');
pk.asa = +(pk.substance =='ASA');
pk.constant = ones(size(pk,1),1);
pk(pk.year>2008 , : ) = [];

demand = NestedLogitDemand(pk);
demand.var.marketsize = 'BL_Unit';
demand.var.price = 'Ptablets_Real'; 
demand.var.nests = 'form substance';
demand.var.quantity = 'Xtablets';
demand.var.market = 'date';
demand.var.panel = 'product';
demand.var.exog = ['marketing1 sw sm  month2 month3 month4 month5 month6 '...
    'month7 month8 month9 month10 month11 month12'];
demand.var.instruments = 'num numg numf numfg numhg numfgh';

demand.settings.paneltype = 'lsdv';
demand.init(); 

results = demand.estimate();
selection = (pk.year==2008 & pk.month==12);
%demand.data = pk(pk.year==2008 & pk.month==12, :);
demand.init();

market = Market(demand);
market.var.firm = 'firm';
market.init();
market.findCosts();

p = demand.data.Ptablets_Real(selection);
c = market.c(selection);

%% Anonymous functions for array manipulation

% Some array operations
times_array = @(x, y) bsxfun(@times, x, y);
plus_array = @(x, y) bsxfun(@plus, x, y);
transpose_array = @(x)  permute(x, [2,1,3]);

% Transform jacobian to array derivative of vector
d_array = @(x) reshape(x, [size(x,1), 1, size(x,2)] );
array2mat = @(x) reshape(x, [size(x,1), size(x,3)] );


%% Share Hessian and FOC Jacobian

% Test code
% A4 = @(x) demand.shares(x) * demand.shares(x)';
% A2 = @(x) bsxfun(@times, 1./(demand.GHGH*demand.shares(x)), ...
%     demand.shares(x) * demand.shares(x)');
% ndA4 = jacobian_array(A4, p);
% ndA2 = jacobian_array(A2, p);

% Depends on:
%       demand.shares(p)
%       demand.shareJacobian(p)
%       demand.alpha
%       demand.sigma
%       demand.GG
%       demand.GHGH
%       market.RR
%       market.c

% Same share jacobian
market.D.initSimulation(168);
demand.initSimulation(168);
s = demand.shares(p);

nds = jacobian(@(x) demand.shares(x), p);
ds = demand.shareJacobian(p);

% Calculate derivative of outer product dA4
dA4 = bsxfun(@times, d_array(ds), s');
dA4 = dA4 + transpose_array(dA4);
 
% Calculate dA2 and dA3 (similar)
sf = (1/(1-demand.sigma(1)) - 1/(1-demand.sigma(2)) );
shInv = 1./(demand.HH*s);
dpart_1 = bsxfun(@times, shInv, dA4);
der = -bsxfun(@times, shInv .^ 2, demand.HH*ds);
dpart_2 = bsxfun(@times, d_array(der), s*s');
dA2 = times_array( sf*demand.HH, dpart_1 + dpart_2);

sf = (1/(1-demand.sigma(2)) - 1 );
shInv = 1./(demand.GG*s);
dpart_1 = bsxfun(@times, shInv, dA4);
der = -bsxfun(@times, shInv .^ 2, demand.GG*ds);
dpart_2 = bsxfun(@times, d_array(der), s*s');
dA3 = times_array( sf*demand.GG, dpart_1 + dpart_2);

%Combine to Hessian
hess = demand.alpha*(-1/(1-demand.sigma(1)) * diag_array(ds) + ...
    dA2 + dA3 + dA4);

% Jacobian of FOC is given by
dfoc = arraymult(times_array(market.RR, hess), (p-c)) + ...
    times(market.RR , ds) + ds;

%% Comparison with numerical calculations

nhess = jacobian_array( @(x) demand.shareJacobian(x), p);

% Small difference with numerical Hessian
max(max(abs(nhess - hess)));

% Numerical jacobian of FOC
ndfoc = jacobian( ...
    @(x) times(market.RR, demand.shareJacobian(x)) * (x - c) + ...
    demand.shares(x), p);

% Small maximum distance between numerical and analytical FOC Jacobian
max(max(abs(ndfoc - dfoc)))
