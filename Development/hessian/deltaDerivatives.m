function [jac, hess] = deltaDerivatives(vx, delta, rc_sigma)
% Calculates both delta Jacobian and Hessian in BLP
% Note that matrix inversion is used rather than solving, due to
% involvement of array mult. As the number of layers is theta, it would be
% best not to solve.
[s, S] = share(vx, delta, rc_sigma);
[dsdd, dsdr] = shareJacobians(vx, delta, rc_sigma);
[d2sdd2,d2sdr2,d2sddr,d2sdrd] = share_hessians(S, vx);

jac = -dsdd\dsdr;
invdsdd = -inv(dsdd);
analder = arraymult(arraymult(d2sdd2, jac, 2), jac) + ...
    arraymult(d2sdrd, jac,2) + arraymult(d2sddr, jac) + d2sdr2;
hess = arraymult(invdsdd, analder);
end

