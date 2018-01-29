function draws = binornd(k, p, n, m)
    draws = rand(n, m) < p;
end