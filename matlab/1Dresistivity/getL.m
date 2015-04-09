function L = getL(n)

    n = n-1;
    h = 1/(n-1);
    L = spdiags(ones(n,1)*[-1 1]/h,[0:1],n-1,n);
   
end