%% setup
n  = 51;
f  = 5;
x  = linspace(0,1,n);
h  = x(2) - x(1);
mt = ones(n,1)./1.^2;

%% get matrix
omega = 2*pi*f;

D  = spdiags(ones(n,1)*[-1 1]/h,[0:1],n-1,n);
w  = ones(n,1); w([1 end]) = .5;
At = omega^2*diags(w.*mt) - D'*D;

[Vr,Dr] = eig(full(At'*At)); [Dr,Ir] = sort(diag(Dr),'descend'); Vr = Vr(:,Ir);
kappar  = Dr(1)/Dr(end);

a  = [1e-1 1 1e1];
L  = [1 10 20];
ls = {'ro','bx','g^'};

P{1} = speye(n);
P{2} = fft(eye(n))/sqrt(n);

for k = 1:1
    for l = 1:length(L)  
        %I = randperm(n); I  = I(1:L(l));
        I = round(linspace(2,n-1,L(l)));
        Pl = P{k}(:,I);

        mu = eigmax(@(x)At'\(Pl*Pl'*(At\x)),n);

        figure;semilogy(1:n,Dr,'k*','markersize',10);hold on;
        xlabel('n');ylabel('eigenvalue');xlim([1 n])

        for i = 1:length(a)
            lambda = a(i)*real(mu);
            B = (At'*At) + (1/lambda)*(Pl*Pl');
            %[Vp,Dp] = eig(full(diags(1./Dr)*ex(Vr'*B*Vr))); [Dp,Ip] = sort(diag(Dp),'descend'); Vp = Vp(:,Ip);
            [Vp,Dp] = eig(full(B)); [Dp,Ip] = sort(diag(Dp),'descend'); Vp = Vp(:,Ip);
            kappap(l,i) = Dp(1)/Dp(end);
            lower = 1/(1+L(l)/(lambda*Dr(end)));
            upper =   (1+L(l)/(lambda*Dr(1)));

            semilogy(1:n,abs(Dp),ls{i},'markersize',10);hold on;
        end
        legend('reduced',['\lambda = 0.1'],['\lambda = 1.0'],['\lambda = 10'],'location','SouthWest');
    end

    kappa{k} = kappap/kappar;
end

savefig(1:3,'../../doc/figs/example3');

latextable(kappa{1},'Horiz',{'$\lambda = 0.1$','$\lambda = 1$','$\lambda = 10$'},'Vert',{['L = ' num2str(L(1))],['L = ' num2str(L(2))],['L = ' num2str(L(3))] },'Hline',[1 NaN],'format','%1.2e','name','../../doc/figs/example3.tex');

