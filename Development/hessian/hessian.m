function hess = hessian( func, x)
% Calculates numerical Hessian of func
dev = diag(.00001*max(abs(x),1e-8*ones(size(x))));
j0 = jacobian(func, x );
hess = zeros(size(j0,1), size(j0,2),length(x));
for i = 1:length(x)
    hess(:,:,i) = (jacobian(func, x + dev(:, i)) - j0) / dev(i, i);
end

end