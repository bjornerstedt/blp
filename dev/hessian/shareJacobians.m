function [dsdelta, dssigma] = shareJacobians( vx, delta, rc_sigma)
[~, si] = share(vx, delta, rc_sigma);
    wsi = si ./ size(si,2); % Row means
dssigma = zeros(size(si, 1), length(rc_sigma));
for k = 1:length(rc_sigma)
    svx = wsi .* vx(:,:,k); % Can be done outside loop
    dssigma(:, k) = sum(svx, 2) - si * sum(svx)';
end
dsdelta = diag(sum(wsi, 2)) - wsi*si';
dsdelta1 = dsdelta(:, 1);
end
   