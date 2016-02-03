a=[2 4  2 4 5; 1 2 2 2 1]';
g=grp2idx(a(:,1))
b=a*[10^3 1]'
grp2idx(b)
dummyvar(g)