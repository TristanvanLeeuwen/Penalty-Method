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

D = spdiags(ones(n,1)*[-1 1]/h,[0:1],n-1,n);
s = [0;ones(n-2,1);0];
A = 1i*2*pi*f*diags(s) + D'*diags(m)*D;



