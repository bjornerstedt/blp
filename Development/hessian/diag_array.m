function d = diag_array(x)
    % Create diagonal array from matrix
    % Should perhaps be sparse array
    [m,n] = size(x);
    d = zeros(m,m,n);
    for i = 1:n;
        d(:,:,i) = diag(x(:, i));
    end;
end
