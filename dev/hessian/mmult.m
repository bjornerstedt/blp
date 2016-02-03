% MMULT multiplies a 3d matrix by a matrix or vector. 
% Without the optional argument, matrix operations are on each matrix 
% defined by the last index (a left multiplication). With anything as
% a third argument, the multiplication is over the first index of 
% either the A or B matrix, depending on which is 3d. 
% Both A and B cannot be 3d as that would imply a 4d resultant matrix.
function C = mmult(A, B, varargin)
rmult = nargin > 2; % Left multiplication by default, with arg right mult
dA = size(A);
dB = size(B);
if size(dB) == 3 & size(dA) == 3
    error('3d by 3d multiplication not implemented')
end
if size(dA,2) == 3
    if rmult
        A = permute(A , [2,3,1])
    end
    if dB(2) == 1 % Unneccessary unless automatic collapse to 2d matrix
        C = zeros([dA(1),1,dA(3)]);
    else
        C = zeros([dA(1), dB(2),dA(3)]);
    end
    for k = 1:dA(3)
        if size(dB, 2) == 2
            C(:,:,k) = A(:,:,k)*B;
        else
            C(:,k) = A(:,:,k)*B;
        end
    end
else
    if size(dB,2) == 3
        if rmult
            B = permute(B , [2,3,1])
        end
        if dA(1) == 1
            C = zeros([1,dB(2),dB(3)]);
        else
            C = zeros([dA(1), dB(2),dB(3)]);
        end
        for k = 1:dB(3)
            if size(dA, 2) == 2
                C(:,:,k) = A*B(:,:,k);
            else
                C(:,k) = A*B(:,:,k);
            end
        end
    else % Both A and B are regular matrices / vectors
        C = A*B;
    end
end
if rmult
    C = permute(C , [3, 1, 2])
end

%{
if lmult
else
    if size(dA,2) == 3
        A = permute(A , [2,3,1])
        if dB(2) == 1
            C = zeros(dA(1),dA(2),1);
        else
            C = zeros([dA(1),dA(2), dB(2)]);
        end
        for k = 1:dA(1)
            if size(dB, 2) == 2
                z = A(k,:,:)
                C(k,:,:) = z*B;
            else
                C(k,:) = A(k,:,:)*B;
            end
        end
    else
        if size(dB,2) == 3
            if dA(1) == 1
                C = zeros(dB(2), dB(3));
            else
                C = zeros([dA(1), dB(2), dB(3)]);
            end
            for k = 1:dB(1)
                if size(dA, 2) == 2
                    C(k,:,:) = A*B(k,:,:);
                else
                    C(k,:) = A*B(k,:,:);
                end
            end
        else % Both A and B are regular matrices / vectors
            C = A*B;
        end
    end
end
%}