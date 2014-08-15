%% setup
f  = [1 3 6 10];
zr = {60};
xr = {100:100:10000};
ns = length(zr{1})*length(xr{1});
nf = length(f);

%% data
vt = dlmread('../../data/marm_20.dat');
mt = 1./vt(:).^2;
n  = size(vt);
h  = [20 20];
X  = (n-1).*h;

Dt = zeros(ns,ns,nf);
Q  = speye(ns);
P  = getP(h,n,zr,xr);    
for k = 1:nf
    At = getA(f(k),mt,h,n);
    Dt(:,:,k) = full(P*(At\(P'*Q)));
end

%% inversion
vref = dlmread('../../data/marm_20.dat'); 
h = [20 20];

z = 0:h(1):X(1);
x = 0:h(2):X(2);
n = [length(z) length(x)];

[zz,xx] = ndgrid(z,x);

v0 = (1.5 + 1*1e-3*max(zz-300,0));
m0 = 1./v0(:).^2;
P  = getP(h,n,zr,xr);

for k = 1:nf
    A0    = getA(f(k),m0,h,n);
    mu(k) = real(eigmax(@(x)A0'\(P'*P*(A0\x)),prod(n)));
end
mask = ones(n); mask(z<200,:) = 0; mask(end,:) = 0;mask(:,[1 end]) = 0;

model.f = f;
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.mref = 1./vref(:).^2;
model.mask = mask(:);

opts.maxit  = 20;
opts.M      = 5;
opts.tol    = 1e-2;
opts.lintol = 1e-1;
opts.method = 'lbfgs';

alpha = 0;

% reduced
mr = m0;
infor = [];
for k = 1:nf
    model.f = f(k);
    fh = @(m)phi(m,Q,Dt(:,:,k),alpha,model);
    [mr,infok] = QGNewton(fh,mr,opts);
    infor = [infor;infok];
end

% penalty
m1 = m0;
info1 = [];
for k = 1:nf
    lambda1 = 1e-1*mu(k);
    model.f = f(k);
    fh = @(m)phi_lambda(m,Q,Dt(:,:,k),alpha,lambda1,model);
    [m1,infok] = QGNewton(fh,m1,opts);
    info1 = [info1;infok];
end

% 
m2 = m0;
info2 = [];
for k = 1:nf
    lambda2 = 1e-0*mu(k);
    model.f = f(k);
    fh = @(m)phi_lambda(m,Q,Dt(:,:,k),alpha,lambda2,model);
    [m2,infok] = QGNewton(fh,m2,opts);
    info2 = [info2;infok];
end

%
m3 = m0;
info3 = [];
for k = 1:nf
    lambda3 = 1e1*mu(k);
    model.f = f(k);
    fh = @(m)phi_lambda(m,Q,Dt(:,:,k),alpha,lambda3,model);
    [m3,infok] = QGNewton(fh,m3,opts);
    info3 = [info3;infok];
end

%% plot
plot2 = @(m)imagesc(x,z,reshape(m,n));

figure;plot2(-P'*P*ones(prod(n),1));colormap(gray);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');

figure;plot2(vref);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');colorbar;caxis([1.5 4.5]);
figure;plot2(v0);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');colorbar;caxis([1.5 4.5]);


figure;
semilogy(sqrt(sum(infor(:,[5,6,7]).^2,2)),'k');hold on;
semilogy(sqrt(sum(info1(:,[5,6,7]).^2,2)),'r');hold on;
semilogy(sqrt(sum(info2(:,[5,6,7]).^2,2)),'b');hold on;
semilogy(sqrt(sum(info3(:,[5,6,7]).^2,2)),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10');
xlabel('iteration');ylabel('||\nabla L||_2');

figure;
plot(infor(:,[8])/infor(1,[8]),'k');hold on;
plot(info1(:,[8])/info1(1,[8]),'r');hold on;
plot(info2(:,[8])/info2(1,[8]),'b');hold on;
plot(info3(:,[8])/info3(1,[8]),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10');
xlabel('iteration');ylabel('||m-m^*||_2');

figure;plot2(real(1./sqrt(mr)));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([1.5 4.5]);
figure;plot2(real(1./sqrt(m1)));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([1.5 4.5]);
figure;plot2(real(1./sqrt(m2)));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([1.5 4.5]);
figure;plot2(real(1./sqrt(m3)));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([1.5 4.5]);

savefig(1:9,'../../doc/figs/2D_exp2');
% 
table = [[1; 2].*infor(end,[1 2])' info1(end,[1 2])' info2(end,[1 2])' info3(end,[1 2])'];
latextable(table,'Horiz',{'reduced','$\lambda = 0.1$','$\lambda = 1$','$\lambda = 10$'},'Vert',{'iterations','PDE solves'},'Hline',[1 NaN],'format','%d','name','../../doc/figs/2D_exp2.tex');
