function [f,g,H,opt] = phi_lambda(m,Q,D,alpha,lambda,model)
% Evaluate penalty misfit
%
% use:
%   [f,g,H] = phi_lambda(m,Q,D,lambda,model)
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
%U = [A;(1/sqrt(lambda))*Q']\[Q;(1/sqrt(lambda))*D];
U = (lambda*(A'*A) + (Q*Q'))\(Q*D + lambda*A'*Q);

%% adjoint field
V = lambda*(A*U - Q);

%% compute f
f = .5*norm(Q'*U - D,'fro')^2 + .5*lambda*norm(A*U - Q,'fro')^2 + .5*alpha*norm(L*m)^2;

%% compute g
g = alpha*(L'*L)*m;

for k = 1:ns
    g = g + real(G(U(:,k))'*V(:,k));
end

%% get H
H = @(x)Hmv(x,m,Q,U,alpha,lambda,model);

%% optimality
opt = [norm(g),  norm(A'*V - Q*(D - Q'*U),'fro'), norm(A*U - Q,'fro'), 0, 0];

end

function y = Hmv(x,m,Q,U,alpha,lambda,model)
%%
n = length(m)+1;
ns = size(Q,2);

%% get matrices
L = getL(n);
A = getA(model.f,m);
G = @(u)getG(m,u);

%% compute mat-vec
y =  alpha*(L'*L)*x;

for k = 1:ns
    y = y + lambda*real(G(U(:,k))'*G(U(:,k))*x - G(U(:,k))'*A*(((1/lambda)*(Q*Q') + (A'*A))\(A'*G(U(:,k))*x)));
end

end
