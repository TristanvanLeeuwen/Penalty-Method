function K = getK(m,u)
% Define Jacobian matrix K(m,u) = d(A(m)'*u)/dm
%
% use:
%   K = getK(f,m,u,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   u - wavefield
%   n - number of gridpoints 
%
% output:
%   G - sparse matrix
%
n = length(m)+1;
h = 1/(n-1);

D = spdiags(ones(n,1)*[-1 1]/h,[0:1],n-1,n);

K = D'*diags(D*u);
