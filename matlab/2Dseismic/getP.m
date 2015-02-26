function P = getP(h,n,zt,xt)
% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zt,xt - arrays defining sampling points
%
% output
%   P     - sparse matrix


z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);

%[zz,xx] = ndgrid(z,x);
% for k = 1:length(zt)
%     i(k) = find((zz(:)==zt(k))&(xx(:)==xt(k)));
% end
% 
% I = speye(prod(n));
% P = I(:,i)/sqrt(prod(h));

P = getLA2(z(:),x(:),[zt(:) xt(:)])'/sqrt(prod(h));

end
function A = getLA2(x1,y1,X2)
    % interior points
    ik = find(X2(:,1) >= x1(1) & X2(:,1) <= x1(end) & X2(:,2) >=y1(1) & X2(:,2) <= y1(end));
    % sizes
    nx1 = length(x1); ny1 = length(y1);
    n1  = nx1*ny1;
    n2  = size(X2,1);
    nk  = length(ik);
    % check
    if ~nk
        A = [];
        return;
    end
    % index sets for sparse matrix
    I = zeros(4*nk,1); J = I; S = J;
    % loop
    l = 1;
    for i = 1:nk
        k = ik(i);
        ix = min(find(x1<=X2(k,1), 1, 'last' ), nx1 - 1);
        iy = min(find(y1<=X2(k,2), 1, 'last' ), ny1 - 1);
        a = ix     + nx1*(iy - 1);
        b = ix     + nx1*(iy);
        c = ix + 1 + nx1*(iy - 1);
        d = ix + 1 + nx1*(iy);
        
        I(l)   = k; J(l)   = a; S(l)   = (X2(k,1) - x1(ix+1))*(X2(k,2) - y1(iy+1))/((x1(ix) - x1(ix+1))*(y1(iy) - y1(iy+1)));
        I(l+1) = k; J(l+1) = b; S(l+1) = (X2(k,1) - x1(ix+1))*(X2(k,2) - y1(iy))  /((x1(ix) - x1(ix+1))*(y1(iy+1) - y1(iy)));
        I(l+2) = k; J(l+2) = c; S(l+2) = (X2(k,1) - x1(ix))  *(X2(k,2) - y1(iy+1))/((x1(ix+1) - x1(ix))*(y1(iy) - y1(iy+1)));
        I(l+3) = k; J(l+3) = d; S(l+3) = (X2(k,1) - x1(ix))  *(X2(k,2) - y1(iy))  /((x1(ix+1) - x1(ix))*(y1(iy+1) - y1(iy)));
        l = l + 4;
    end
    % construct matrix
    A = sparse(I(1:l-1),J(1:l-1),S(1:l-1),n2,n1);
end
