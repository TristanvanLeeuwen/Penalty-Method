function [f,g,H] = phi_lambda(m,Q,D,lambda,model)
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
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% forward solve
U = [sqrt(lambda)*A;P]\[sqrt(lambda)*P'*Q;D];

%% compute f
f = .5*norm(P*U - D,'fro')^2 + .5*lambda*norm(A*U - P'*Q,'fro')^2;

%% compute g
g = zeros(size(m));

for k = 1:ns
    g = g + real(lambda*G(U(:,k))'*(A*U(:,k) - P'*Q(:,k)));
end

%% get H
H = @(x)Hmv(x,m,Q,U,lambda,model);

end

function y = Hmv(x,m,Q,U,lambda,model)
%%
ns = size(Q,2);

%% get matrices
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% compute mat-vec
y = zeros(size(m));

for k = 1:ns
    y = y + real(lambda*G(U(:,k))'*G(U(:,k))*x + lambda^2*G(U(:,k))'*A*((P'*P + lambda*(A'*A))\(A'*G(U(:,k))*x)));
end

end
