function [mu,v] = eigmax(fh,n)

    delta = 1;
    tol   = 1e-6;
    iter  = 0;
    
    mu = 0;    
    v  = randn(n,1);

    while (delta>tol)&&(iter<100)
        vt  = fh(v);
        mut = vt'*v;

        delta = abs(mut - mu)/abs(mut);
        mu    = mut;
        v     = vt/norm(vt);
        
        iter  = iter + 1;
    end
end