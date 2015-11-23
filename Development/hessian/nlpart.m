function expmu = nlpart(vx, rc_sigma)
% Computes the mu or nonlinear part of utility
[n,m,K] = size(vx);
mu = zeros(n, m);
for k = 1:K
    mu = mu + vx(:,:,k) * rc_sigma(k);
end
expmu = exp(mu);
end