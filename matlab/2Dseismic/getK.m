function K = getK(f,m,v,h,n)
% Define Jacobian matrix K(m,u) = d(A(m)'*v)/dm
%
% use:
%   K = getK(f,m,u,d,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   u - wavefield
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   K - sparse matrix
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

K = omega^2*spdiags(1e-6*w.*v,0,N,N) - 1i*omega*spdiags(1e-3*(1-w).*v./(2*sqrt(m)),0,N,N);

