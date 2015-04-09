%% setup
f    = 5;
mfun = @(x)1 + exp(-1e1*(x(:)-.5).^2);

alpha = 1e-6;
sigma = 0;
delta = linspace(-1,1,50);
%% data
n  = 201;
x  = linspace(0,1,n-1);
mt = mfun(x);

At = getA(f,mt);
Qt = 1e1*speye(n)*sqrt(n-1);
Qt = Qt(:,[1 end]);
Dt = Qt'*(At\(Qt));

noise = randn(size(Dt)); noise = sigma*norm(Dt(:))*noise/norm(noise(:));
SNR   = 20*log10(norm(Dt(:))/norm(noise(:)))
Dt    = Dt + noise;

%% inversion
n  = 101;
x  = linspace(0,1,n-1);
Q  = 1e1*speye(n)*sqrt(n-1);
Q  = Q(:,[1 end]);
model.f = f;

mt = mfun(x);
m0 = ones(n-1,1);
A0 = getA(f,m0);
mu = real(eigmax(@(x)A0'\(Q*Q'*(A0\x)),n));

opts.maxit  = 100;
opts.M      = 100;
opts.tol    = 1e-9;
opts.lintol = 1e-3;
opts.method = 'GN';

opts.isreal = 1;


% reduced
[fr0,gr,Hr] = phi(mt,Q,Dt,alpha,model);
[Vr,Dr,flagr] = eigs(Hr,n-1,6,'LM',opts); [Dr,Ir] = sort(diag(Dr),'descend'); Vr = Vr(:,Ir);

for k = 1:length(delta)
    for l = 1:length(delta)
        fr(k,l) = phi(mt + delta(k)*Vr(:,1) + delta(l)*Vr(:,2),Q,Dt,alpha,model);
        %qr(k,l) = fr0 + delta(k)*(gr'*Vr(:,1)) + delta(l)*(gr'*Vr(:,2)) + .5*delta(k)^2*Dr(1) + .5*delta(l)^2*Dr(2);
    end
end


