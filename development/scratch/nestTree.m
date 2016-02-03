function nt  = nestTree( types )
    nt = (1:types(end))';
    for i = (length(types)-1):-1:1
        nt = [ reshape(repmat(1:types(i), size(nt, 1), 1), [], 1), ...
            repmat(nt, types(i), 1)];
    end
end

