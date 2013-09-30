%%
% test for the introduction of auxiliary variables

% min_{x} 1/2*|y-A*x|^2 + |B(x)|_1       
% min_{z = (x,u)} 1/2*|y-A*x|^2 + |u|_1 + i_{B(x)=u}

addpath('../');
addpath('../toolbox/');

n = 100;
p = 50;
q = 2*n;

A = randn(p,n);
y = randn(p,1);
B = randn(q,n);


% |x-x'| + gamma* |y-A*x'|^2
% x'-x + A'(A*x'-y) = 0

% projector on {(x,u) \ u=B(x)}
if 0
    ProjTmp = @(x,u,tu)deal(x-B'*tu,u+tu);
    Proj  = @(x,u)ProjTmp(x,u, (eye(size(B,1),size(B,1))+B*B')\(B*x-u) );
else
    ProjTmp = @(x,tx)deal(tx,B*tx);
    Proj = @(x,u)ProjTmp(x, (eye(size(B,2),size(B,2))+B'*B)\(B'*u+x) );
end


concat = @(x,u)[x(:); u(:)];
X = @(z)z(1:n);
U = @(z)z(n+1:n+q);

ProxF = @(x,gamma)(eye(n)+gamma*A'*A)\(x+gamma*A'*y);
ProxG = @(u,gamma)perform_soft_thresholding(u, gamma);
% prox of constraint
MyProj = @(z)Proj(X(z),U(z));
ProxA  = @(z,gamma)apply_multiple_ouput( concat, MyProj,z );
% prox of initial functional
ProxB  = @(z,gamma)concat( ProxF(X(z),gamma), ProxG(U(z),gamma) );



options.niter = 5000;
options.report = @(z)1/2*norm(y-A*X(z))^2 + norm(B*X(z),1);
[x,R] = perform_dr(zeros(n+q,1), ProxB, ProxA, options);


clf; 
plot(log10(R(1:end/2) - min(R))); axis tight;

