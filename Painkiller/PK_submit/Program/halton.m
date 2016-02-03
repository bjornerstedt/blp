function f=halton(n,s)
k=round(log(n+1)./log(s));
phi=0;
i=1;
while i<=k
      x=phi;
      j=1;
      while j<s
            y=phi+(j/s^i);
            x=[x; y];
            j=j+1;
      end
      phi=x;
      i=i+1;
end

x=phi;
j=1;
while j<s && size(x,1)<(n+1)
    y=phi+(j/s^i);
    x=[x; y];
    j=j+1;
end
f=x(2:n+1,1);