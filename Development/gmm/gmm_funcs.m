function beta = gmm_funcs(demand)
% Various GMM estimates
%display('Simple Estimate')
%beta = opt1(demand)
    alpha = [-.3]';

display('Nested Estimate')
beta = opt2(demand, alpha)

display('Simultaneous Estimate')

market = Market(demand);
market.var.firm = 'productid';
market.settings.paneltype = 'none';
market.var.exog = 'w';
market.init();
    
beta = market.gmm_estimate( alpha)

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

function alpha = opt2(demand, alpha)
    % Nested finding of alpha
    W = inv(demand.Z' * demand.Z);
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

