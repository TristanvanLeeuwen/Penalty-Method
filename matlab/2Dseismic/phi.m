function [f,g,H,opt] = phi(m,D,alpha,model)
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
if isfield(model,'mref')
    mref = model.mref;
else
    mref = 0*m;
end
if isfield(model,'mask')
    mask = model.mask;
else
    mask = 1;
end
%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
Q = getP(model.h,model.n,model.zs,model.xs);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% forward solve
U = A\Q;

%% compute f
f = .5*norm(P'*U - D,'fro')^2 + .5*alpha*norm(L*m)^2;

%% adjoint solve
V = A'\(P*(D - P'*U));

%% compute g
g = alpha*(L'*L)*m;

for k = 1:size(U,2)
    g = g + real(G(U(:,k))'*V(:,k));
end
g = mask.*g;

%% get H
H = @(x)Hmv(x,m,U,alpha,model);

%% optimality
opt = [norm(g),  norm(A'*V - P*(D - P'*U),'fro'), norm(A*U - Q,'fro'), norm(m-mref), norm((D - P'*U),'fro')];

end

function y = Hmv(x,m,U,alpha,model)
%%
if isfield(model,'mask')
    mask = model.mask;
else
    mask = 1;
end
%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% compute mat-vec
y = mask.*x;
y = alpha*(L'*L)*y;

for k = 1:size(U,2);
    y = y + real(G(U(:,k))'*(A'\((P*P')*(A\(G(U(:,k))*x)))));
end
y = mask.*y;
end
