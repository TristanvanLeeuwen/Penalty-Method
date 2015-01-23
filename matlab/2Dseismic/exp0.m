%% setup
f  = 10;
c0 = 2;
c1 = 2.5;
c2 = 2.25;


zr = [[20:20:980]']';
xr = [20*ones(49,1)]';
zs = [[20:20:980]']';
xs = [980*ones(49,1)]';

mfun = @(zz,xx)1./(c0 + (c1-c0)*exp(-5e-5*(xx(:)-300).^2 - 5e-5*(zz(:)-300).^2) + (c2-c0)*exp(-5e-5*(xx(:)-700).^2 - 5e-5*(zz(:)-700).^2)).^2;

sigma = 0;
alpha = 0;
delta = 10*linspace(-1,1,20);

%% data
n  = [51 51];
h  = [20 20];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

mt = mfun(zz,xx);
At = getA(f,mfun(zz,xx),h,n);
Pt = getP(h,n,zr,xr);
Qt = getP(h,n,zs,xs);
Dt = Pt'*(At\Qt);

noise = randn(size(Dt)); noise = sigma*noise*(norm(Dt(:))/norm(noise(:)));
SNR   = 20*log10(norm(Dt(:))/norm(noise(:)))
Dt    = Dt + noise;

%% inversion
n  = [51 51]; N = prod(n);
h  = [20 20];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

mask = ones(n); mask([1:1 end:end],:) = 0;mask(:,[1 end:end]) = 0;

mt = mfun(zz,xx);

model.f = f;
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.zs = zs;
model.xs = xs;
model.mref = mt;
model.mask = mask(:);

m0 = ones(prod(n),1)/c0.^2;
A0 = getA(f,m0,h,n);
P  = getP(h,n,zr,xr);
mu = real(eigmax(@(x)A0'\(P*P'*(A0\x)),prod(n)));

opts.isreal = 1;

%v1 = sin(pi*zz(:)/1000); v1 = v1/norm(v1);
%v2 = sin(2*pi*zz(:)/1000);  v2 = v2/norm(v2);
v1 = 1./(c0 + (c1-c0)*exp(-5e-5*(xx(:)-300).^2 - 5e-5*(zz(:)-300).^2)).^2 - m0;
v2 = 1./(c0 + (c2-c0)*exp(-5e-5*(xx(:)-700).^2 - 5e-5*(zz(:)-700).^2)).^2 - m0;

% reduced
[fr0,gr,Hr] = phi(m0,Dt,alpha,model);
%[Vr,Dr,flagr] = eigs(Hr,N,6,'LM',opts); [Dr,Ir] = sort(diag(Dr),'descend'); Vr = Vr(:,Ir);
Vr = [v1 v2];
for k = 1:length(delta)
    for l = 1:length(delta)
        fr(k,l) = phi(m0 + delta(k)*Vr(:,1) + delta(l)*Vr(:,2),Dt,alpha,model);
        %qr(k,l) = fr0 + delta(k)*(gr'*Vr(:,1)) + delta(l)*(gr'*Vr(:,2)) + .5*delta(k)^2*Dr(1) + .5*delta(l)^2*Dr(2);
    end
end

% penalty
lambda1 = 1e-1*mu;
lambda2 = 1e-0*mu;
lambda3 = 1e1*mu;

[fp10,gp1,Hp1] = phi_lambda(m0,Dt,alpha,lambda1,model);  
%[Vp1,Dp1,flag1] = eigs(Hp1,N,6,'LM',opts); [Dp1,Ip1] = sort(diag(Dp1),'descend'); Vp1 = Vp1(:,Ip1);
Vp1 = [v1 v2];
for k = 1:length(delta)
    for l = 1:length(delta)
        fp1(k,l) = phi_lambda(m0 + delta(k)*Vp1(:,1) + delta(l)*Vp1(:,2),Dt,alpha,lambda1,model);
        %qp1(k,l) = fp10 + delta(k)*(gp1'*Vp1(:,1)) + delta(l)*(gp1'*Vp1(:,2)) + .5*delta(k)^2*Dp1(1) + .5*delta(l)^2*Dp1(2);
    end
end

[fp20,gp2,Hp2] = phi_lambda(m0,Dt,alpha,lambda2,model);
%[Vp2,Dp2,flag2]     = eigs(Hp2,N,6,'LM',opts); [Dp2,Ip2] = sort(diag(Dp2),'descend'); Vp2 = Vp2(:,Ip2);
Vp2 = [v1 v2];
for k = 1:length(delta)
    for l = 1:length(delta)
        fp2(k,l) = phi_lambda(m0 + delta(k)*Vp2(:,1) + delta(l)*Vp2(:,2),Dt,alpha,lambda2,model);
        %qp2(k,l) = fp20 + delta(k)*(gp2'*Vp2(:,1)) + delta(l)*(gp2'*Vp2(:,2)) + .5*delta(k)^2*Dp2(1) + .5*delta(l)^2*Dp2(2);
    end
end

lambda = 1e1*mu;
[fp30,gp3,Hp3] = phi_lambda(m0,Dt,alpha,lambda3,model);
%[Vp3,Dp3,flag3]     = eigs(Hp3,N,6,'LM',opts); [Dp3,Ip3] = sort(diag(Dp3),'descend'); Vp3 = Vp3(:,Ip3);
Vp3 = [v1 v2];
for k = 1:length(delta)
    for l = 1:length(delta)
        fp3(k,l) = phi_lambda(m0 + delta(k)*Vp3(:,1) + delta(l)*Vp3(:,2),Dt,alpha,lambda3,model);
        %qp3(k,l) = fp30 + delta(k)*(gp3'*Vp3(:,1)) + delta(l)*(gp3'*Vp3(:,2)) + .5*delta(k)^2*Dp3(1) + .5*delta(l)^2*Dp3(2);
    end
end

save('exp0');
%% plot
load('exp0')
plot2 = @(m)imagesc(x,z,reshape(m,n));

figure;plot2(Vr(:,1));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.05 .05])
figure;plot2(Vp1(:,1));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.05 .05])
figure;plot2(Vp2(:,1));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.05 .05])
figure;plot2(Vp3(:,1));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.05 .05])
% 
figure;plot2(Vr(:,2));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.07 .07])
figure;plot2(Vp1(:,2));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.07 .07])
figure;plot2(Vp2(:,2));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.07 .07])
figure;plot2(Vp3(:,2));axis equal tight;ylabel('x_1 [m]');xlabel('x_2 [m]');caxis([-.07 .07])

figure;contourf(delta,delta,fr,20);xlabel('\delta_2');ylabel('\delta_1');
figure;contourf(delta,delta,fp1,20);xlabel('\delta_2');ylabel('\delta_1');
figure;contourf(delta,delta,fp2,20);xlabel('\delta_2');ylabel('\delta_1');
figure;contourf(delta,delta,fp3,20);xlabel('\delta_2');ylabel('\delta_1');


%savefig(1:12,'../../doc/figs/2D_exp0');
