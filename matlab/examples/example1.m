A = @(m)m;
G = @(m,u)u;
R = 1;
P = 1;
d = 0;
q = 1;

misfit = @(m,u).5*norm(P'*u-d)^2 + .5*norm(R*m)^2;

g = @(w)[G(w(1),w(2))'*w(3) + R'*R*w(1); A(w(1))'*w(3) + P*(P'*w(2) - d); A(w(1))*w(2) - q];
H = @(w)[R'*R, G(w(1),w(3)), G(w(1),w(2))'; G(w(1),w(3))' P*P' A(w(1))'; G(w(1),w(2)) A(w(1)) 0];

% g = @(w)[w(2)*w(3) + w(1); w(2)-2 + w(1)*w(3); w(1)*w(2)-1];
% H = @(w)[w(1) w(3) w(2);w(3) 1 w(1);w(2) w(1) 0];

w0(1) = .5;
w0(2) = A(w0(1))\q;
w0(3) = A(w0(1))\(P*(d - P'*w0(2)));
w0 = w0(:);

% all-at-once
w1  = w0;
ng1 = norm(g(w1));
k   = 1;
while (ng1(end)>1e-10)&&(k<100)
    gk     = g(w1(:,k)); 
    Hk     = H(w1(:,k));
    w1(:,k+1) = w1(:,k)-Hk\gk;   
    k = k + 1;
    ng1(k) = norm(g(w1(:,k)));
end

%% reduced
w2  = w0;
ng2 = norm(g(w2));
k   = 1;
while (ng2(end)>1e-10)&&(k<50)
    gk      = g(w2(:,k)); 
    Hk      = H(w2(:,k));
    gred    = gk(1);
    Hred    = Hk(1,1) - Hk(1,[2:3])*inv(Hk([2:3],[2:3]))*Hk([2:3],1);
    w2(1,k+1) = w2(1,k) - Hred\gred;
    w2(2,k+1) = A(w2(1,k+1))\q;
    w2(3,k+1) = A(w2(1,k+1))\(P*(d - P'*w2(2,k+1)));
    k = k + 1;
    ng2(k)  = norm(g(w2(:,k)));
    f2(k)   = misfit(w2(1,k),w2(2,k));
end

%% penalty
lambda = 1e3;
w3  = w0;
ng3 = norm(g(w3));
w3(3,1) = lambda*(A(w3(1,1))*w3(2,1) - q);
k   = 1;
while (ng3(end)>1e-10)&&(k<200)
    gk = [G(w3(1,k),w3(2,k))'*w3(3,k) + R'*R*w3(1,k);P*(P'*w3(1,k) - d) + A(w3(1,k))'*w3(3,k)];
    Hk = [lambda*G(w3(1,k),w3(2,k))'*G(w3(1,k),w3(2,k)) + R'*R ,lambda*G(w3(1,k),w3(2,k))'*A(w3(1,k)) + G(w3(1,k),w3(3,k)) ;lambda*A(w3(1,k))'*G(w3(1,k),w3(2,k)) + G(w3(1,k),w3(3,k)) ,P*P' + lambda*A(w3(1,k))'*A(w3(1,k))];
    w3([1 2],k+1) = w3([1 2],k) - Hk\gk;
    w3(3,k+1) = lambda*(A(w3(1,k+1))*w3(2,k+1) - q);
    k = k + 1;
    ng3(k)  = norm(g(w3(:,k)));
     f3(k)  = misfit(w3(1,k),w3(2,k)) + lambda*norm(A(w3(1,k))*w3(2,k) - q)^2;
end


%% plot
ut = linspace(0,10,100);
mt = linspace(0,2,100);

for k = 1:length(ut)
    for l = 1:length(mt)
        fh(k,l) = misfit(mt(l),ut(k));
        ph(k,l) = .5*lambda*norm(A(mt(l))*ut(k) - q)^2;
    end
end

for l = 1:length(mt)
    u(l) = A(mt(l))\q;
end

figure;
contourf(mt,ut,fh,50);hold on;plot(mt,u,'w--');
hold on;plot(w1(1,:),w1(2,:),'w*-',w2(1,:),w2(2,:),'w^-');

figure;
contourf(mt,ut,fh+ph,50);hold on;plot(mt,u,'w--');hold on;
plot(w3(1,:),w3(2,:),'wo-');

figure;semilogy(1:length(ng1),ng1,1:length(ng2),ng2,1:length(ng3),ng3)