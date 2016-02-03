% Calculates the share Hessian arrays as function of delta, 
% arrays of dimensions: MxMxM, MxKxM,MxMxK and MxKxK where
% M is the number of products in S and K the number of of random coefs.
% NOTE: The function works, but is slow as it loops over individuals
% Can thus benefit from C++ calculation. Using symmetry in the outer
% product also reduces the number of multiplication to a half in 2d and
% more in 3d.
function [hd, hr, hdr, hrd] = share_hessians(S, vx)
M = size(S, 1);
V = permute(vx, [1,3,2]);

for i = 1:size(S, 2)
    hess = zeros(M, M, M);
    s = S(:, i);
    ss = s*s';
    own = s - 2 * s .* s; % Alt: own=s with no inner if statement
    for j = 1:M
        hess(:,:,j) = 2 * ss * s(j);
        for k = 1:M
            hess(k, k, j) = hess(k, k, j) - ss(j,k); % deltaJac*s(j);
            if k ~= j
                hess(j, k, j) = hess(j, k, j) - ss(j,k);
                hess(k, j, j) = hess(k, j, j) - ss(j,k); 
            end
        end
        hess(j,j,j) = hess(j,j,j) + own(j);
    end
    hessdr = arraymult(hess, V(:,:,i), 2);
    hessrd = arraymult(hess, V(:,:,i));
    if i == 1
        hd = hess; % quad goes here 
        hdr = hessdr; % quad goes here 
        hrd = hessrd; % quad goes here 
        hr = arraymult(hessdr, V(:,:,i));
    else
        hd = hd + hess;
        hdr = hdr + hessdr; % quad goes here 
        hrd = hrd + hessrd; % quad goes here 
        hr = hr + arraymult(hessdr, V(:,:,i));
    end
end
hd = hd/size(S, 2); % Quadrature requires weighting
hr = hr/size(S, 2);
hdr = hdr/size(S, 2);
hrd = hrd/size(S, 2);