function newedel = finddelta(s, vx, rc_sigma, delta0)
i = 0;
maxdev = 100;

tolerance = 1e-14;
fpmaxit = 1000;
edelta = exp(delta0); 
expmu = nlpart(vx, rc_sigma); % Creates obj.mu

while maxdev > tolerance && i < fpmaxit
    eg = bsxfun(@times, edelta, expmu );
    denom = 1 ./ ( 1 + sum(eg));

        share = mean(bsxfun(@times, eg, denom), 2); % Row means

    newedel = edelta .* s ./ share;
    maxdev = max(abs(newedel - edelta));
    %                     dev = newedel - edeltasel;   % Is this motivated?
    %                     maxdev = dev'*dev;
    edelta = newedel;
    i = i + 1;
end
end
