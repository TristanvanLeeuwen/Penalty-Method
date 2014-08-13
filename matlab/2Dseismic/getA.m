function A = getA(f,m,h,n)
% Define 2D Helmholtz matrix with absorbing boundary conditions
%
% use:
%   A = getA(f,m,h,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   A - sparse matrix
%

%%
omega = 2*pi*f;
N     = prod(n);

%% Stiffness matricex
D1 = spdiags(ones(n(1),1)*[1 -2 1]/h(1).^2,[-1:1],n(1),n(1)); D1(1,1:2) = [1 -1]/h(1); D1(end,end-1:end) = [-1 1]/h(1);
D2 = spdiags(ones(n(2),1)*[1 -2 1]/h(2).^2,[-1:1],n(2),n(2)); D2(1,1:2) = [1 -1]/h(2); D2(end,end-1:end) = [-1 1]/h(2);
S  = kron(speye(n(2)),D1) + kron(D2,speye(n(1)));

% D1 = spdiags(ones(n(1),1)*[-1 1]/h(1),[0:1],n(1)-1,n(1));
% D2 = spdiags(ones(n(2),1)*[-1 1]/h(2),[0:1],n(2)-1,n(2));
% S  = kron(speye(n(2)),-D1'*D1) + kron(-D2'*D2,speye(n(1)));

%% Mass matrix, including ABC
a = ones(n); a(:,[1 end]) = 0; a([1 end],:) = 0; a = a(:);
%b = zeros(n); b(:,[1 end]) = 1./h(1); b([1 end],:) = 1./h(2); b = b(:);
b = 1-a;
M = omega^2*spdiags(1e-6*a.*m,0,N,N) + 1i*omega*spdiags(b.*sqrt(1e-6*m),0,N,N);

%% 
A = M + S;



