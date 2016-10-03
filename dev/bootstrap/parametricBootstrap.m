% bootstrap returns an array with bootreps draws of beta.
function bootpars = parametricBootstrap(beta, varcovar, bootreps)
    % Numerical errors cause varcovar to be nonsymmetric:
    varcovar = (varcovar + varcovar')/2;
    bootpars = mvnrnd(beta, varcovar, bootreps);
end
