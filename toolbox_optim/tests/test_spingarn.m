%%
% Test for multiple function using Spingarn digonal trick.
% Simple test problem:
%   min_x i_{A*x=y}(x) + |B(x)|_1 + |C(x)|_1
% where B and C are orthogonal matrices.

addpath('../');
addpath('../toolbox/');


n = 100;
p = 40;
A = randn(p,n);
y = randn(p,1);
[B,R] = qr(randn(n));
[C,R] = qr(randn(n));

FA = @(x)norm(A*x-y);
FB = @(x)norm(B*x,1);
FC = @(x)norm(C*x,1);
options.report = @(x)struct('Const', FA(x), 'Obj', FB(x)+FC(x));

pA = A'*(A*A')^(-1);
ProxFA = @(x,tau)x + pA*(y-A*x);
ProxFB = @(x,tau)B'*perform_soft_thresholding(B*x, tau);
ProxFC = @(x,tau)C'*perform_soft_thresholding(C*x, tau);

options.niter = 5000;
[x,R] = perform_dr_spingarn(zeros(n,1), {ProxFA ProxFB ProxFC}, options);

Const = s2v(R, 'Const');
Obj = s2v(R, 'Const');

clf; 
subplot(2,1,1);
plot(log10(Const/Const(1))); axis tight;
subplot(2,1,2);
plot(log10(Obj(1:end/2) - min(Obj))); axis tight;