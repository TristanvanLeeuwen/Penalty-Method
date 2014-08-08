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
ns = size(Q,2);

%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% forward solve
%U = [sqrt(lambda)*A;P]\[sqrt(lambda)*P'*Q;D];
U = (lambda*(A'*A) + (P'*P))\(P'*D + lambda*A'*P'*Q);

%% adjoint field
V = lambda*(A*U - P'*Q);

%% compute f
f = .5*norm(P*U - D,'fro')^2 + .5*lambda*norm(A*U - P'*Q,'fro')^2 + .5*alpha*norm(L*m)^2;

%% compute g
g = alpha*(L'*L)*m;

for k = 1:ns
    g = g + real(G(U(:,k))'*V(:,k));
end

%% get H
H = @(x)Hmv(x,m,Q,U,alpha,lambda,model);

%% optimality
opt = [norm(g),  norm(A'*V - P'*(D - P*U),'fro'), norm(A*U - P'*Q,'fro')];

end

function y = Hmv(x,m,Q,U,alpha,lambda,model)
%%
ns = size(Q,2);

%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% compute mat-vec
y =  alpha*(L'*L)*x;

for k = 1:ns
    y = y + real(lambda*G(U(:,k))'*G(U(:,k))*x - lambda^2*G(U(:,k))'*A*((P'*P + lambda*(A'*A))\(A'*G(U(:,k))*x)));
end

end
