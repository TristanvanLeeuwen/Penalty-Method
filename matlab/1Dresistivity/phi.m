function [f,g,H,opt] = phi(m,Q,D,alpha,model)
% Evaluate reduced misfit
%
% use:
%   [f,g,H] = phi(m,Q,D,model)
%
% input:
%
% output:
%
%

%%
n = length(m)+1;
ns = size(Q,2);

%% get matrices
L = getL(n);
A = getA(model.f,m);
G = @(u)getG(m,u);

%% forward solve
U = A\(Q);

%% compute f
f = .5*norm(Q'*U - D,'fro')^2 + .5*alpha*norm(L*m)^2;

%% adjoint solve
V = A'\(Q*(D - Q'*U));

%% compute g
g = alpha*(L'*L)*m;

for k = 1:ns
    g = g + real(G(U(:,k))'*V(:,k));
end

%% get H
H = @(x)Hmv(x,m,Q,U,alpha,model);

%% optimality
opt = [norm(g),  norm(A'*V - Q*(D - Q'*U),'fro'), norm(A*U - Q,'fro'), 0, 0];

end

function y = Hmv(x,m,Q,U,alpha,model)
%%
n = length(m)+1;
ns = size(Q,2);

%% get matrices
L = getL(n);
A = getA(model.f,m);
G = @(u)getG(m,u);

%% compute mat-vec
y = alpha*(L'*L)*x;

for k = 1:ns
    y = y + real(G(U(:,k))'*(A'\((Q*Q')*(A\(G(U(:,k))*x)))));
end

end
