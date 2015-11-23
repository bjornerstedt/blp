function hess = crosshessian( func, x, y)
% Calculates the numerical cross derivative of f(x, y)
dev = diag(.00001*max(abs(y),1e-8*ones(size(y))));
deriv = @(w) jacobian(@(z) func(z, w), x );
j0 = deriv( y );
hess = zeros(size(j0,1), length(x),length(y));
for i = 1:length(y)
    hess(:,:,i) = (deriv(y + dev(:, i)) - j0) / dev(i, i);
end

end