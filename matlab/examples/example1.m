%% problem setup
n = 1;
alpha = 1;
A = @(m)1/(27*m.^3 - 54*m.^2 + 36.1*m - 7.9);
G = @(m,u)-u*(3*27&m^2 + 2*54*m + 36.1)/(27*m.^3 - 54*m.^2 + 36.1*m - 7.9)^2;
K = @(m,v)-v*(3*27&m^2 + 2*54*m + 36.1)/(27*m.^3 - 54*m.^2 + 36.1*m - 7.9)^2;
R = @(m,u,v);
P = 1;
d = 1;
q = 1;

% n = 2;
% alpha = .1;
% A = @(m)[1 + 2*m(1), .5; .5, 1 + m(2)];
% G = @(m,u)diags([2;1].*u);
% K = @(m,v)diags([2;1].*v);
% R = @(m,u,v)0*speye(n);
% P = [1 0;0 1];
% q = [1;1];
% d = P*(A([1;1])\q);

%% gradient and Hessian
g = @(w)[G(w(1:n),w(n+1:2*n))'*w(2*n+1:3*n) + alpha*w(1:n); A(w(1:n))'*w(2*n+1:3*n) + P*(P'*w(n+1:2*n) - d); A(w(1:n))*w(n+1:2*n) - q];
H = @(w)[R(w(1:n),w(n+1:2*n),w(2*n+1:3*n)) + alpha*speye(n), K(w(1:n),w(2*n+1:3*n))',G(w(1:n),w(n+1:2*n))';K(w(1:n),w(2*n+1:3*n)), P*P',A(w(1:n))'; G(w(1:n),w(n+1:2*n)), A(w(1:n)), 0*speye(n)];

gred = @(w)G(w(1:n),w(n+1:2*n))'*w(2*n+1:3*n) + alpha*w(1:n);
Hred = @(w)G(w(1:n),w(n+1:2*n))'*inv(A(w(1:n))')*(P*P')*inv(A(w(1:n)))*G(w(1:n),w(n+1:2*n)) + alpha*speye(n) - K(w(1:n),w(2*n+1:3*n))'*inv(A(w(1:n)))*G(w(1:n),w(n+1:2*n)) - G(w(1:n),w(n+1:2*n))*inv(A(w(1:n))')*K(w(1:n),w(2*n+1:3*n)) + R(w(1:n),w(n+1:2*n),w(2*n+1:3*n));

%% initial guess, feasible point
m0 = [.5];
u0 = A(m0)\q;
v0 = A(m0)'\(P*(d - P'*u0));

%% all-at-once
k   = 1;
w1  = [m0;u0;v0]; 
ng1 = norm(g(w1));
while (ng1(end)>1e-10)&&(k<100)
    dw = -H(w1(:,k))\g(w1(:,k));
    
    w1(:,k+1) = w1(:,k) + dw;
    
    ng1(k+1)  = norm(g(w1(:,k+1)));
    k = k + 1;
end

%% reduced
k   = 1;
w2  = [m0;u0;v0]; 
ng2 = norm(g(w2));
while (ng2(end)>1e-10)&&(k<100)
    dm = -Hred(w2(:,k))\gred(w2(:,k));
    
    w2(0*n+1:1*n,k+1) = w2(1:n,k) + dm;
    w2(1*n+1:2*n,k+1) = A(w2(1:n,k+1))\q;
    w2(2*n+1:3*n,k+1) = A(w2(1:n,k+1))'\(P*(d - P'*w2(n+1:2*n,k+1)));
    
    ng2(k+1)        = norm(g(w2(:,k+1)));
    k = k + 1;
end 

%% plot

% at  = linspace(0,1,50);
% phi = zeros(length(at));
% for k = 1:length(at)
%     for l = 1:length(at)
%         mt = at([k l])';
%         ut = A(mt)\q;
%         phi(k,l) = .5*norm(P*ut - d)^2 + alpha*0.5*norm(mt)^2;
%     end
% end

ut  = linspace(0,2,50);
mt  = linspace(0,2,50);
misfit = zeros(length(ut),length(mt));
phi    = zeros(length(mt));
constraint = zeros(length(mt));
for l = 1:length(mt)
    ul = A(mt(l))\q;
    constraint(l) = ul;
    phi(l)        = .5*norm(P*ul - d)^2 + alpha*0.5*norm(mt(l))^2;
    for k = 1:length(ut)
        misfit(k,l) = .5*norm(P*ut(k) - d)^2 + alpha*0.5*norm(mt(l))^2;
    end
end

figure;
semilogy(1:length(ng1),ng1,1:length(ng2),ng2)

figure;
contourf(mt,ut,misfit,50);hold on;
plot(mt,constraint,'k--');hold on;
plot(w1(1,:),w1(2,:),'wo-',w2(1,:),w2(2,:),'kx-');

figure;
plot(mt,phi)