function P = getP(h,n,zt,xt)
% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zt,xt - cell arrays defining receiver arrays, e.g., zt{1} = .1; xt{1} = [.1:.1:.9]';
%
% output
%   P     - sparse matrix

I1 = speye(n(1));
I2 = speye(n(2));
z  = h(1)*[0:n(1)-1]';
x  = h(2)*[0:n(2)-1]';

P  = kron(I2(ismember(x,xt{1}),:),I1(ismember(z,zt{1}),:));

for k = 2:length(zt)
    P  = [P;kron(I2(ismember(x,xt{k}),:),I1(ismember(z,zt{k}),:))];
end
