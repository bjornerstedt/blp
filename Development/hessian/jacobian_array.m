function J = jacobian_array(func, x0)
    % Numerical jacobian, extension for matrix valued function
    % From: http://fabcol.free.fr/pdf/lectnotes4.pdf
    %
    f = feval(func, x0);
    k = length(x0);
    [m,n] = size(f);
    J = zeros(m,n,k);
    dev = diag(.00001*max(abs(x0),1e-8*ones(size(x0))));
    for i = 1:k;
        ff = feval(func, x0+dev(:,i));
        J(:,:,i) = (ff-f)/dev(i,i);
    end;
end
