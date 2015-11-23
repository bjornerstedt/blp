function J = jacobian(func, x0)
    % Numerical jacobian, optionally used in MixedLogitDemand
    % From: http://fabcol.free.fr/pdf/lectnotes4.pdf
    % [df]=numgrad(func,x0,method,param)
    %
    x0 = x0(:);
    f = feval(func,x0);
    m = length(x0);
    n = length(f);
    J = zeros(n,m);
    dev = diag(.00001*max(abs(x0),1e-8*ones(size(x0))));
    for i = 1:m;
        ff = feval(func, x0+dev(:,i));
        J(:,i) = (ff-f)/dev(i,i);
    end;
end

%{
    % method = 'c' -> centered difference
    % = 'l' -> left difference
    %
    elseif (lower(method)=='c')
        for i=1:m;
            ff= feval(func,x0+dev(:,i));
            fb= feval(func,x0-dev(:,i));
            J(:,i) = (ff-fb)/(2*dev(i,i));
        end;
    else
        error('Bad method specified')
    end
%}