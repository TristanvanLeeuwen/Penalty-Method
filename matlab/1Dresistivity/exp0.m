%% setup
n  = 11;
f  = 5;
x  = linspace(0,1,n-1);
mt = 1 + exp(-1e1*(x-.5).^2);
%mt = ones(n-1,1);mt(1:floor(n/2))=2;

%% data
At = getA(f,mt);
Q  = speye(n);
Q  = Q(:,[2 end-1]);
Dt = Q'*(At\(Q));

M = size(Q,2);
%% factorization
m0 = ones(n-1,1);
A0 = getA(f,m0);
mu = eigmax(@(x)A0'\(Q*Q'*(A0\x)),n);

lambda = 1e-1*real(mu);


[V,D]   = eig(full(A0'*A0)); [D,I] = sort(diag(D)); V = V(:,I);
[Vt,Dt] = eig(full((A0'*A0) + (1/lambda)*(Q*Q'))); [Dt,It] = sort(diag(Dt)); Vt = Vt(:,It);

figure; plot(1:n,D,1:n,Dt)

figure; plot(1:n,Dt-D,1:n,ones(n,1)/lambda);

sum((Dt-D)*lambda)

[(1/(1+M/lambda/D(1)))*D(end)/D(1)  Dt(end)/Dt(1) (1+M/lambda/D(end))*D(end)/D(1)] 

for k = 2:n
    a = sign(V(:,k)'*Vt(:,k));
    subplot(2,5,k-1);plot(1:n,V(:,k),1:n,a*Vt(:,k));
end