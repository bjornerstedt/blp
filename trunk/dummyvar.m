function res = dummyvar(vec)
    vec = vec - min(vec) + 1;
    res = zeros(length(vec), max(vec) );
    for i = 1:length(vec)
        res(i, vec(i) ) = 1;
    end
end