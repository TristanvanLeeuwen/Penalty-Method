function R = getR(f,m,u,v,h,n)
% Define Jacobian matrix R(m,u,v) = d(G(m,u)*v)/dm
%
% use:
%   R = getR(f,m,u,v,d,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   u - wavefield
%   v - adjoint wavefield
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   G - sparse matrix
%

%%
omega = 2*pi*f;
N     = prod(n);

%%
w = [0 ones(1,n(1)-2) 0];
if n(2)>1
    w = w(:)*[0 ones(1,n(2)-2) 0];
end
w = w(:);

R = .25*1i*omega*spdiags(1e-3*(1-w).*conj(u).*v.*(m.^(-3/2)),0,N,N);

