%% setup
n  = [51 51];
N  = prod(n);
h  = [20 20];
f  = 5;
c0 = 2;
c1 = 2.5;
c2 = 2.25;
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);
zr = {[60:80:940],[60 940]};
xr = {[60 940],[140:80:860]};
ns = length(zr{1})*length(xr{1}) + length(zr{2})*length(xr{2});

vt = c0 + (c1-c0)*exp(-5e-5*(xx-300).^2 - 5e-5*(zz-300).^2) + (c2-c0)*exp(-5e-5*(xx-700).^2 - 5e-5*(zz-700).^2);
mt = 1./vt(:).^2;

%% data
At = getA(f,mt,h,n);
P  = getP(h,n,zr,xr);
Q  = speye(ns);
Dt = P*(At\(P'*Q));

%% inversion
model.f = f;
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;

m0 = ones(prod(n),1)/c0.^2;
A0 = getA(f,m0,h,n);
mu = real(eigmax(@(x)A0'\(P'*P*(A0\x)),prod(n)));

opts.maxit  = 100;
opts.M      = 5;
opts.tol    = 1e-9;
opts.lintol = 1e-1;
opts.method = 'GN';

alpha = 1e1;

% reduced
fh = @(m)phi(m,Q,Dt,alpha,model);
[mr,infor] = QGNewton(fh,m0,opts);

% penalty
lambda = 1e-1*mu;
fh = @(m)phi_lambda(m,Q,Dt,alpha,lambda,model);
[m1,info1] = QGNewton(fh,m0,opts);

lambda = 1e-0*mu;
fh = @(m)phi_lambda(m,Q,Dt,alpha,lambda,model);
[m2,info2] = QGNewton(fh,m0,opts);

lambda = 1e1*mu;
fh = @(m)phi_lambda(m,Q,Dt,alpha,lambda,model);
[m3,info3] = QGNewton(fh,m0,opts);

%% plot
% figure;
% semilogy(infor(:,1),infor(:,[5]),'k--',infor(:,1),infor(:,[7]),'k-');legend('L_m','L_v');hold on;
% semilogy(info1(:,1),info1(:,[5]),'r--',info1(:,1),info1(:,[7]),'r-');hold on
% semilogy(info2(:,1),info2(:,[5]),'b--',info2(:,1),info2(:,[7]),'b-');hold on;
% semilogy(info3(:,1),info3(:,[5]),'g--',info3(:,1),info3(:,[7]),'g-');

plot2 = @(m)imagesc(x,z,reshape(m,n));

figure;plot2(-P'*P*ones(prod(n),1));colormap(gray);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');

figure;plot2(mt);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');colorbar;

figure;
semilogy(sqrt(sum(infor(:,[5,6,7]).^2,2)),'k');hold on;
semilogy(sqrt(sum(info1(:,[5,6,7]).^2,2)),'r');hold on;
semilogy(sqrt(sum(info2(:,[5,6,7]).^2,2)),'b');hold on;
semilogy(sqrt(sum(info3(:,[5,6,7]).^2,2)),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','southwest');
xlabel('iteration');ylabel('||\nabla L||_2');


figure;plot2(mt);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(mr);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m1);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m2);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m3);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');

savefig(1:8,'../../doc/figs/2D_exp1');

table = [[1; 2].*infor(end,[1 2])' info1(end,[1 2])' info2(end,[1 2])' info3(end,[1 2])'];
latextable(table,'Horiz',{'reduced','$\lambda = 0.1$','$\lambda = 1$','$\lambda = 10$'},'Vert',{'iterations','PDE solves'},'Hline',[1 NaN],'format','%d','name','../../doc/figs/2D_exp1.tex');
