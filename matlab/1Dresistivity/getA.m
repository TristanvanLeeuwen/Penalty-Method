function A = getA(f,m)
% Define system matrix
%
% use:
%   A = getA(f,m,h,n)
%
% input:
%   f - frequency [Hz]
%   m - parameter
%
% output:
%   A - sparse matrix
%
n = length(m)+1;
h = 1/(n-1);
omega = 2*pi*f;

D = spdiags(ones(n,1)*[-1 1]/h,[0:1],n-1,n);
w = ones(n,1); w([1 end]) = 0.5;
A = 1i*omega*diags(w) + D'*diags(m)*D;



