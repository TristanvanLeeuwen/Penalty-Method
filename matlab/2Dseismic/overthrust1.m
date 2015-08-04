%% setup
v = dlmread('overthrust_25.dat');
v = 1e-3*v(1:2:end,1:2:end);
dx = 50;
f  = 2;

xr = [[100:200:19900]']';
zr = [100*ones(100,1)]';
xs = [[200:200:19800]']';
zs = [100*ones(99,1)]';

sigma = 0;
alpha = 5;
maxit = 50;

%% data
n  = size(v);
h  = [dx dx];
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
v = dlmread('overthrust_25.dat');
v = 1e-3*v(1:4:end,1:4:end);
dx = 100;

n  = size(v);
h  = [dx dx];
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

m0 = 1./v.^2;
for k = 1:20
   m0=filter2([0.25 .5 0.25]'*[0.25 .5 0.25],padarray(m0,[1,1],'replicate'),'valid');
end
m0 = m0(:);

A0 = getA(f,m0,h,n);
P  = getP(h,n,zr,xr);
Q  = getP(h,n,zs,xs);
D0 = P'*(A0\Q);
mu = real(eigmax(@(x)A0'\(P*P'*(A0\x)),prod(n)));

opts.maxit  = maxit;
opts.M      = 10;
opts.tol    = 1e-6;
opts.lintol = 1e-1;
opts.method = 'lbfgs';

% reduced
fh = @(m)phi(m,Dt,alpha,model);
[mr,infor] = QGNewton(fh,m0,opts);
Ar = getA(f,mr,h,n);
Dr = P'*(Ar\Q);

% penalty
lambda = 1e-2*mu;
fh = @(m)phi_lambda(m,Dt,alpha,lambda,model);
[m1,info1] = QGNewton(fh,m0,opts);
A1 = getA(f,mr,h,n);
D1 = P'*(A1\Q);

lambda = 1e-0*mu;
fh = @(m)phi_lambda(m,Dt,alpha,lambda,model);
[m2,info2] = QGNewton(fh,m0,opts);
A2 = getA(f,mr,h,n);
D2 = P'*(A2\Q);

lambda = 1e2*mu;
fh = @(m)phi_lambda(m,Dt,alpha,lambda,model);
[m3,info3] = QGNewton(fh,m0,opts);
A3 = getA(f,mr,h,n);
D3 = P'*(A3\Q);

save('overthrust1');

%% plot
load('overthrust1');

plot2 = @(m)imagesc(1e-3*x,1e-3*z,reshape(m,n),[0.02 .18]);

figure;plot2(mref);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');colorbar;hold on; plot(1e-3*xs(5:10:end),1e-3*zs(5:10:end),'k*',1e-3*xr(1:10:end),1e-3*zr(1:10:end),'kv','markersize',10,'linewidth',2)

figure;
semilogy(sqrt(sum(infor(:,[5,6,7]).^2,2)),'k');hold on;
semilogy(sqrt(sum(info1(:,[5,6,7]).^2,2)),'r');hold on;
semilogy(sqrt(sum(info2(:,[5,6,7]).^2,2)),'b');hold on;
semilogy(sqrt(sum(info3(:,[5,6,7]).^2,2)),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
xlabel('iteration');ylabel('||\nabla L||_2');axis square tight

figure;
plot(infor(:,8),'k');hold on;
plot(info1(:,8),'r');hold on;
plot(info2(:,8),'b');hold on;
plot(info3(:,8),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
xlabel('iteration');ylabel('||m^k - m^*||_2');axis square tight

figure;
semilogy(infor(:,9),'k');hold on;
semilogy(info1(:,9),'r');hold on;
semilogy(info2(:,9),'b');hold on;
semilogy(info3(:,9),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
xlabel('iteration');ylabel('||d - P^Tu^k||_2');axis square tight

figure;
plot(infor(:,9),infor(:,7),'k-o');hold on;
plot(info1(:,9),info1(:,7),'r-o');hold on;
plot(info2(:,9),info2(:,7),'b-o');hold on;
plot(info3(:,9),info3(:,7),'g-o');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10','location','northeast');
ylabel('||A(m)u^k - d||_2');xlabel('||P^Tu^k - d||_2');axis tight;ylim([-1e-4 4e-3]);xlim([0 3]);axis square

figure;plot2(m0);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(mr);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m1);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m2);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');
figure;plot2(m3);axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');

Is = 50;
Ir = 1:2:length(xr);
figure;plot(1e-3*xr(Ir),real(Dt(Ir,Is)),'k--',1e-3*xr,real(D0(:,Is)),'k-');ylim([-.1 .1]);set(gca,'plotboxaspectratio',[3 1 1]);xlabel('x_r [km]');ylabel('\Re(d)');
figure;plot(1e-3*xr(Ir),real(Dt(Ir,Is)),'k--',1e-3*xr,real(Dr(:,Is)),'k-');ylim([-.1 .1]);set(gca,'plotboxaspectratio',[3 1 1]);xlabel('x_r [km]');ylabel('\Re(d)');
figure;plot(1e-3*xr(Ir),real(Dt(Ir,Is)),'k--',1e-3*xr,real(D1(:,Is)),'k-');ylim([-.1 .1]);set(gca,'plotboxaspectratio',[3 1 1]);xlabel('x_r [km]');ylabel('\Re(d)');
figure;plot(1e-3*xr(Ir),real(Dt(Ir,Is)),'k--',1e-3*xr,real(D2(:,Is)),'k-');ylim([-.1 .1]);set(gca,'plotboxaspectratio',[3 1 1]);xlabel('x_r [km]');ylabel('\Re(d)');
figure;plot(1e-3*xr(Ir),real(Dt(Ir,Is)),'k--',1e-3*xr,real(D3(:,Is)),'k-');ylim([-.1 .1]);set(gca,'plotboxaspectratio',[3 1 1]);xlabel('x_r [km]');ylabel('\Re(d)');

savefig(1:15,'../../doc/figs/2D_overthrust1');

table = [[1; 2].*infor(end,[1 2])' info1(end,[1 2])' info2(end,[1 2])' info3(end,[1 2])'];
latextable(table,'Horiz',{'reduced','$\lambda = 0.1$','$\lambda = 1$','$\lambda = 10$'},'Vert',{'iterations','PDE solves'},'Hline',[1 NaN],'format','%d','name','../../doc/figs/2D_overthrust1.tex');
