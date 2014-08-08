%% setup
n  = [51 51];
N  = prod(n);
h  = [20 20];
f  = 7;
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
v0=randn(N,1);for k = 1:100, v1 = A0'\(P'*P*(A0\v0));mu=real(v1'*v0);v0=v1/norm(v1);end

alpha = 1e1;
delta = 10*linspace(-1,1,20);
opts.isreal = 1;

% reduced
[fr,gt,Hr] = phi(mt,Q,Dt,alpha,model);
[Vr,Dr,flagr]    = eigs(Hr,N,6,0,opts);
for k = 1:length(delta)
    fr(k) = phi(mt + delta(k)*Vr(:,1),Q,Dt,alpha,model);
end

% penalty
lambda1 = 1e-1*mu;
lambda2 = 1e-0*mu;
lambda3 = 1e1*mu;

[fp1,gp1,Hp1] = phi_lambda(mt,Q,Dt,alpha,lambda1,model);
[Vp1,Dp1,flag1]     = eigs(Hp1,N,6,0,opts);
for k = 1:length(delta)
    fp1(k) = phi_lambda(mt + delta(k)*Vp1(:,1),Q,Dt,alpha,lambda1,model);
end

[fp2,gp2,Hp2] = phi_lambda(mt,Q,Dt,alpha,lambda2,model);
[Vp2,Dp2,flag2]     = eigs(Hp2,N,6,0,opts);
for k = 1:length(delta)
    fp2(k) = phi_lambda(mt + delta(k)*Vp2(:,1),Q,Dt,alpha,lambda2,model);
end

lambda = 1e1*mu;
[fp3,gp3,Hp3] = phi_lambda(mt,Q,Dt,alpha,lambda3,model);
[Vp3,Dp3,flag3]     = eigs(Hp3,N,6,0,opts);
for k = 1:length(delta)
    fp3(k) = phi_lambda(mt + delta(k)*Vp3(:,1),Q,Dt,alpha,lambda3,model);
end

%% plot
figure;
for k = 1:6
    subplot(2,3,k);imagesc(reshape(real(Vr(:,k)),n));title(num2str(Dr(k,k)));axis equal tight
end

figure;
for k = 1:6
    subplot(2,3,k);imagesc(reshape(real(Vp1(:,k)),n));title(num2str(Dp1(k,k)));axis equal tight
end

figure;
for k = 1:6
    subplot(2,3,k);imagesc(reshape(real(Vp2(:,k)),n));title(num2str(Dp2(k,k)));axis equal tight
end

figure;
for k = 1:6
    subplot(2,3,k);imagesc(reshape(real(Vp3(:,k)),n));title(num2str(Dp3(k,k)));axis equal tight
end

figure;
plot(delta,fr,'k--',delta,fp1,'r',delta,fp2,'b',delta,fp3,'g');
