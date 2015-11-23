% vx matrix: (product x ind x rc)
function [f, indshare] = share(vx, delta, rc_sigma)
edelta = exp(delta); 
eg = bsxfun(@times, edelta, nlpart(vx, rc_sigma));
denom = 1 ./ ( 1 + sum(eg));
indshare = bsxfun(@times, eg , denom);
if false
    f = indshare * quadw; % Row means
else
    f = mean(indshare, 2); % Row means
end
end

