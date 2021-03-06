%% setup
v = dlmread('../../data/bp_50.dat');
v = 1e-3*v(1:2:201,200 + [1:2:401]);

f  = 1;

xr = [[100:100:19900]']';
zr = [100*ones(199,1)]';
xs = [[100:100:19900]']';
zs = [100*ones(199,1)]';

sigma = 0;
alpha = 0.1;

%% data
n  = [101 201];
h  = [100 100];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

At = getA(f,1./v(:).^2,h,n);
Pt = getP(h,n,zr,xr);
Qt = getP(h,n,zs,xs);
Dt = Pt'*(At\Qt);

noise = randn(size(Dt)); noise = sigma*noise*(norm(Dt(:))/norm(noise(:)));
SNR   = 20*log10(norm(Dt(:))/norm(noise(:)))
Dt    = Dt + noise;

%% inversion
n  = [101 201];
h  = [100 100];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

mref = 1./v(:).^2;

model.f = f;
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.zs = zs;
model.xs = xs;
model.mref = mref;

m0 = 1./(v(1)+min(.6e-7*(zz(:)).^2,1.75+.1*1e-3*zz(:))).^2;
A0 = getA(f,m0,h,n);
P  = getP(h,n,zr,xr);
mu = real(eigmax(@(x)A0'\(P*P'*(A0\x)),prod(n)));

opts.maxit  = 20;
opts.M      = 10;
opts.tol    = 1e-6;
opts.lintol = 1e-1;
opts.method = 'lbfgs';

% reduced
fh = @(m)phi(m,Dt,alpha,model);
[mr,infor] = QGNewton(fh,m0,opts);

% penalty
lambda = 1e-2*mu;
fh = @(m)phi_lambda(m,Dt,alpha,lambda,model);
[m1,info1] = QGNewton(fh,m0,opts);

lambda = 1e-0*mu;
fh = @(m)phi_lambda(m,Dt,alpha,lambda,model);
[m2,info2] = QGNewton(fh,m0,opts);

lambda = 1e2*mu;
fh = @(m)phi_lambda(m,Dt,alpha,lambda,model);
[m3,info3] = QGNewton(fh,m0,opts);

save('bp1');

%% plot
load('bp1');

plot2 = @(m)imagesc(x,z,reshape(1./sqrt(m),n),[1.5 5]);

figure;plot2(mref);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');colorbar;hold on; plot(xs,zs,'k*',xr,zr,'kv','markersize',10,'linewidth',2)

figure;
semilogy(sqrt(sum(infor(:,[5,6,7]).^2,2)),'k');hold on;
semilogy(sqrt(sum(info1(:,[5,6,7]).^2,2)),'r');hold on;
semilogy(sqrt(sum(info2(:,[5,6,7]).^2,2)),'b');hold on;
semilogy(sqrt(sum(info3(:,[5,6,7]).^2,2)),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
xlabel('iteration');ylabel('||\nabla L||_2');axis tight

figure;
plot(infor(:,8),'k');hold on;
plot(info1(:,8),'r');hold on;
plot(info2(:,8),'b');hold on;
plot(info3(:,8),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
xlabel('iteration');ylabel('||m^k - m^*||_2');axis tight

figure;
semilogy(infor(:,9),'k');hold on;
semilogy(info1(:,9),'r');hold on;
semilogy(info2(:,9),'b');hold on;
semilogy(info3(:,9),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
xlabel('iteration');ylabel('||d - P^Tu^k||_2');axis tight

figure;plot2(mr);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m1);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m2);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m3);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');

savefig(1:8,'../../doc/figs/2D_bp1');

table = [[1; 2].*infor(end,[1 2])' info1(end,[1 2])' info2(end,[1 2])' info3(end,[1 2])'];
latextable(table,'Horiz',{'reduced','$\lambda = 0.1$','$\lambda = 1$','$\lambda = 10$'},'Vert',{'iterations','PDE solves'},'Hline',[1 NaN],'format','%d','name','../../doc/figs/2D_exp4.tex');
