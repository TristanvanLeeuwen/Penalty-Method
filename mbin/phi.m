function [f,g,H] = phi(m,Q,D,model)
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
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% forward solve
U = A\(P'*Q);

%% compute f
f = .5*norm(P*U - D,'fro')^2;

%% compute g
V = A'\(P'*(P*U - D));
g = zeros(size(m));

for k = 1:ns
    g = g - real(G(U(:,k))'*V(:,k));
end

%% get H
H = @(x)Hmv(x,m,Q,U,model);

end

function y = Hmv(x,m,Q,U,model)
%%
ns = size(Q,2);

%% get matrices
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% compute mat-vec
y = zeros(size(m));

for k = 1:ns
    y = y + real(G(U(:,k))'*(A'\((P'*P)*(A\(G(U(:,k))*x)))));
end

end
