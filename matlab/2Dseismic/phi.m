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
ns = size(Q,2);

%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% forward solve
U = A\(P'*Q);

%% compute f
f = .5*norm(P*U - D,'fro')^2 + .5*alpha*norm(L*m)^2;

%% adjoint solve
V = A'\(P'*(D - P*U));

%% compute g
g = alpha*(L'*L)*m;

for k = 1:ns
    g = g + real(G(U(:,k))'*V(:,k));
end

%% get H
H = @(x)Hmv(x,m,Q,U,alpha,model);

%% optimality
opt = [norm(g),  norm(A'*V - P'*(D - P*U),'fro'), norm(A*U - P'*Q,'fro')];

end

function y = Hmv(x,m,Q,U,alpha,model)
%%
ns = size(Q,2);

%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% compute mat-vec
y = alpha*(L'*L)*x;

for k = 1:ns
    y = y + real(G(U(:,k))'*(A'\((P'*P)*(A\(G(U(:,k))*x)))));
end

end
