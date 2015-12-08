function beta = gmm_funcs(demand)
% Various GMM estimates
display('Simple Estimate')
beta = opt1(demand)
display('Nested Estimate')
beta = opt2(demand)
display('Simultaneous Estimate')
beta = opt3(demand)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = opt1(demand)
% Simple regression using fminunc
W = inv(demand.Z' * demand.Z);
beta = [.3,.3,.3]';
options = optimoptions(@fminunc, 'Algorithm','quasi-newton', 'MaxIter',50);

[beta,fval,exitflag] = fminunc(@(x)objective(x, demand), beta, options);

    function val = objective(beta, demand)
        xi = demand.y - demand.X*beta;
        xiZ = xi'*demand.Z;
        val = xiZ*W*xiZ';
    end
end

function alpha = opt2(demand)
% Nested finding of alpha
W = inv(demand.Z' * demand.Z);
alpha = [-.3]';
options = optimoptions(@fminunc, 'Algorithm','quasi-newton', 'MaxIter',50);

[alpha,fval,exitflag] = fminunc(@(x)objective(x, demand), alpha, options);

    function val = objective(alpha, demand)
        X0 = demand.X(: , 2:end);
        beta = (X0' * X0)\ (X0' * (demand.y - demand.X(:, 1) * alpha));
        xi = demand.y - demand.X * [alpha; beta];
        xiZ = xi' * demand.Z;
        val = xiZ * W * xiZ';
    end
end

function alpha = opt3(demand)
% Simultaneous estimate of demand and costs over alpha
W = inv(demand.Z' * demand.Z);
market = Market(demand);
market.var.firm = 'productid';
market.settings.paneltype = 'none';
market.var.exog = 'constant';
market.init();
alpha = [-.3]';
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'MaxIter',50);

[alpha,fval,exitflag] = fminunc(@(x)objective(x, demand, market), alpha, options);

    function val = objective(alpha, demand, market)
        % Residuals as function of demand params:
        X0 = demand.X(: , 2:end);
        beta = (X0' * X0)\ (X0' * (demand.y - demand.X(:, 1) * alpha));
        xi = demand.y - demand.X * [alpha; beta];
        
        market.D.alpha = alpha;
        market.findCosts();
        
        % Residuals of market estimation
        gamma = mean(market.c);
        eta = market.c - gamma;
        xiZ = xi' * demand.Z;
        val = xiZ * W * xiZ' + eta' * eta ./ size(market.c, 1);
    end
end
