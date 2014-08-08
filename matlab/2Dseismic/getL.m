function L = getL(h,n)

    D1 = spdiags(ones(n(1),1)*[-1 1]/h(1),[0:1],n(1)-1,n(1));
    D2 = spdiags(ones(n(2),1)*[-1 1]/h(2),[0:1],n(2)-1,n(2));
    
    L  = [kron(speye(n(2)),D1); kron(D2,speye(n(1)))];
end