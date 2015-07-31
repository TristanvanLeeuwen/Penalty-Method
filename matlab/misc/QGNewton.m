function [x,info] = QGNewton(fh,x,options)
% Simple Quasi/Gauss-Newton method with Wolfe linesearch
%
% use:
%   [xn,info] = QGNewton(fh,x0,options)
%
% input:
%   fh - function handle to misfit of the form [f,g] = fh(x)
%        where f is the function value, g is the gradient of the same size
%        as the input vector x. 
%   x0 - initial guess
%
%   options.itermax - max iterations [default 10]
%   options.tol     - tolerance on 2-norm of gradient [1e-3]
%   options.M       - history size [5]
%   options.fid     - file id for output [1]
%   options.write   - save iterates to disk [0]
%   options.method  - {Quasi,Gauss}
%
% output:
%   xn - final estimate
%

if nargin<3
    options = [];
end

% parse parameters
tol     = getoption(options,'tol',1e-6);
lintol  = getoption(options,'lintol',1e-3);
itermax = getoption(options,'maxit',10);
M       = getoption(options,'M',5);
fid     = getoption(options,'fid',1);
write   = getoption(options,'write',0);
method  = getoption(options,'method','lbfgs');

% initialization
n         = length(x);
converged = 0;
iter      = 0;
S         = zeros(n,0);
Y         = zeros(n,0);
lambda    = 1;

% initial evaluation
[f,g,H,opt]  = fh(x);
nfeval = 1;
info   = [iter,nfeval,lambda,f,opt];
fprintf(fid,'# iter, # eval, stepsize, f(x)       , ||L_m||_2,    ||L_u||_2,     ||L_v||_2,   ||m-m_r||_2,  ||Pu - d||_2\n');
formatstr = '%6d, %6d, %1.2e, %1.5e, %1.5e, %1.5e, %1.5e, %1.5e, %1.5e\n';
fprintf(fid,formatstr,info);
if write
    dlmwrite(['x_' num2str(iter) '.dat'],x);
end
% main loop
while ~converged    
    % compute search direction
    switch method
        case 'lbfgs'
            if (iter>0)&&(M>0)
                s = B(-g,S,Y);
                lambda = 1;
            else
                s = -g/norm(g);
            end
        otherwise
            [s,~,~,ncg] = pcg(H,-g,lintol,100);
            nfeval = nfeval + ncg;
    end
    p = -(s'*g)/(g'*g);
    
    if (p < 0)
        fprintf(fid,'Loss of descent, reset history\n');
        S = zeros(n,0);
        Y = zeros(n,0);
        s = B(-g,S,Y);
    end
    
    % linesearch
    [ft,gt,Ht,opt,lambda,lsiter] = wWolfeLS(fh,x,f,g,s,lambda);
    nfeval = nfeval + lsiter;
    
    % update
    if lambda>0
        xt = x + lambda*s;

        S = [S xt - x];
        Y = [Y gt - g];

        if size(S,2)>M
            S = S(:,end-M+1:end);
            Y = Y(:,end-M+1:end);
        end
        f = ft;
        g = gt;
        H = Ht;
        x = xt;
    end
    iter = iter + 1;
    info   = [info;[iter,nfeval,lambda,f,opt]];
    fprintf(fid,formatstr,info(end,:));
    if write
        dlmwrite(['x_' num2str(iter) '.dat'],x);
    end
    
    % check convergence
    converged = (iter>=itermax)||(norm(g)<tol)||(lambda<1e-10);
    
end

end

function z = B(x,S,Y)
% apply lbfgs inverse Hessian to vector
%
% use:
%   z = B(x,S,Y,b0,Binit)
%
% input:
%   x - vector of length n
%   S - history of steps in n x M matrix
%   Y - history of gradient differences in n x M matrix
%
% output
%   z - vector of length n
%

M = size(S,2);

alpha = zeros(M,1);
rho   = zeros(M,1);
for k = 1:M
    rho(k) = 1/(Y(:,k)'*S(:,k));
end
q = x;
% first recursion
for k = M:-1:1
    alpha(k) = rho(k)*S(:,k)'*q;
    q        = q - alpha(k)*Y(:,k);
end

% apply `initial' Hessian
if M>0 
    a = (Y(:,end)'*S(:,end)/(Y(:,end)'*Y(:,end)));
else
    a = 1;
end
z = a*q;
% second recursion
for k = 1:M
    beta = rho(k)*Y(:,k)'*z;
    z    = z + (alpha(k) - beta)*S(:,k);
end
end

function [ft,gt,Ht,opt,lambda,lsiter] = wWolfeLS(fh,x0,f0,g0,s0,lambda)
% Simple Wolfe linesearch, adapted from
% (http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3).
%
%

lsiter = 0;
c1 = 1e-2;
c2 = 0.9;
done = 0;
mu = 0;
nu = inf;
lambda = .5*lambda;
while ~done
    if nu < inf
        lambda = (nu + mu)/2;
    else
        lambda = 2*lambda;
    end
    
    if lsiter < 10
        [ft,gt,Ht,opt] = fh(x0 + lambda*s0);
        lsiter = lsiter + 1;
    else
        lambda = 0;
        break;
    end
    
    %fprintf(1,'      >%d, %1.5e, %1.5e, %1.5e\n',lsiter, lambda, ft, gt'*s0);
    
    if ft > f0 + c1*lambda*g0'*s0
        nu = lambda;
    elseif gt'*s0 < c2*g0'*s0
        mu = lambda;
    else
        done = 1;
    end
end
end
