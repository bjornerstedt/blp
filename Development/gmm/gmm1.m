function beta = gmm1(demand)
W = inv(demand.Z' * demand.Z);
beta = [.3,.3,.3]';
options = optimoptions(@fminunc, ...
    'Algorithm','quasi-newton' ,...
    'MaxIter',50);

[beta,fval,exitflag] = fminunc(@(x)objective1(x, demand), beta, options);

function val = objective1(beta, demand)
xi = demand.y - demand.X*beta;
xiZ = xi'*demand.Z;
val = xiZ*W*xiZ';
end

end