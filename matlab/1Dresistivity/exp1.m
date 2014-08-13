%% setup
n  = 51;
f  = 5;
x  = linspace(0,1,n-1);
mt = 1 + exp(-1e1*(x-.5).^2);
%mt = ones(n-1,1);mt(1:floor(n/2))=2;

%% data
At = getA(f,mt);
Q  = speye(n)*(n-1);
Q  = Q(:,[2 end-1]);
Dt = Q'*(At\(Q));

%% inversion
model.f = f;

m0 = ones(n-1,1);
A0 = getA(f,m0);
mu = real(eigmax(@(x)A0'\(Q*Q'*(A0\x)),n));

opts.maxit  = 100;
opts.M      = 100;
opts.tol    = 1e-9;
opts.lintol = 1e-3;
opts.method = 'GN';

alpha = 1e-6;

% reduced
fh = @(m)phi(m,Q,Dt,alpha,model);
[mr,infor] = QGNewton(fh,m0,opts);

% penalty
lambda1 = 1e-1*mu;
fh = @(m)phi_lambda(m,Q,Dt,alpha,lambda1,model);
[m1,info1] = QGNewton(fh,m0,opts);

lambda2 = 1e0*mu;
fh = @(m)phi_lambda(m,Q,Dt,alpha,lambda2,model);
[m2,info2] = QGNewton(fh,m0,opts);

lambda3 = 1e1*mu;
fh = @(m)phi_lambda(m,Q,Dt,alpha,lambda3,model);
[m3,info3] = QGNewton(fh,m0,opts);


%% plot
% figure;
% semilogy(infor(:,1),infor(:,[5]),'k--',infor(:,1),infor(:,[7]),'k-');legend('L_m','L_v');hold on;
% semilogy(info1(:,1),info1(:,[5]),'r--',info1(:,1),info1(:,[7]),'r-');hold on
% semilogy(info2(:,1),info2(:,[5]),'b--',info2(:,1),info2(:,[7]),'b-');hold on;
% semilogy(info3(:,1),info3(:,[5]),'g--',info3(:,1),info3(:,[7]),'g-');
% xlabel('iteration');ylabel('magnitude of gradient');

figure;
semilogy(sqrt(sum(infor(:,[5,6,7]).^2,2)),'k');hold on;
semilogy(sqrt(sum(info1(:,[5,6,7]).^2,2)),'r');hold on;
semilogy(sqrt(sum(info2(:,[5,6,7]).^2,2)),'b');hold on;
semilogy(sqrt(sum(info3(:,[5,6,7]).^2,2)),'g');
legend('reduced','\lambda = 0.1','\lambda = 1','\lambda = 10');
xlabel('iteration');ylabel('||\nabla L||_2');

figure;
plot(x,mt,'k--',x,mr,'k',x,m1,'r',x,m2,'b',x,m3,'g');
legend('true','reduced','\lambda = 0.1','\lambda = 1','\lambda = 10');
xlabel('x');ylabel('m(x)');

savefig(1:2,'../../doc/figs/1D_exp1');

table = [[1; 2].*infor(end,[1 2])' info1(end,[1 2])' info2(end,[1 2])' info3(end,[1 2])'];
latextable(table,'Horiz',{'reduced','$\lambda = 0.1$','$\lambda = 1$','$\lambda = 10$'},'Vert',{'iterations','PDE solves'},'Hline',[1 NaN],'format','%d','name','../../doc/figs/1D_exp1.tex');
