function z = arraymult( x, y, varargin )
% ARRAYMULT multiply 3d array with matrix/vector
% Multiplication with vector reduces dimensionality
% Only two types of post-multiplication have been defined
% Compare with mmult.m
isarray = @(a)(length(size(a)) == 3);
dotmult = nargin == 3 && varargin{1} == 2;
if isarray(x)
    if dotmult
        x = permute(x, [1,3,2]);
    end
    sz = size(x);
    sz(2) = size(y, 2);
    z = zeros(sz);
    for i = 1:size(x,3)
        z(:, :, i) = x(:, :, i) * y;
    end
elseif isarray(y)
    sz = size(y);
    sz(1) = size(x, 1);
    z = zeros(sz);
    for i = 1:size(y,3)
        z(:, :, i) = x * y(:, :, i);
    end
else
    z = x * y;
end
if dotmult
    z = permute(z, [1,3,2]);
end
if isarray(z)
    if sz(1) == 1
        z = reshape(z, sz([2,3]) );
    elseif sz(2) == 1
        z = reshape(z, sz([1,3]) );
    elseif sz(3) == 1
        z = reshape(z, sz([1,2]) );
    end
end
end

