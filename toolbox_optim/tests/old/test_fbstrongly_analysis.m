%%
% test for resolution of analysis prior denoising
% min_x 1/2*|x-y|^2 + lambda*|K*x|_1
% min_x F(x) + G(K*x)   with   G=lambda*|.|_1,  F=1/2*|.-y|^2

addpath('../');
addpath('../toolbox/');

n = 100;
p = 200;
lambda = 1;
K = randn(p,n);
KS = K';

ProxG = @(u,tau)perform_soft_thresholding(u, lambda*tau);
ProxGS = compute_dual_prox(ProxG);
GradFS = @(x)x-y;
F = @(x)1/2*norm(x-y)^2;
G = @(u)lambda*norm(u,1);
options.FS = F;
options.GS = @(u)norm(u,'inf')/lambda;

L = norm(K)^2;
options.method = 'fista';
options.niter = 2000;
[x, fs, Constraint, Hu] = perform_fb_strongly(zeros(n,1), K, KS, GradFS, ProxGS, L, options);

%% 
% Compute the primal error.

E = [];
for i=1:size(Hu,2)
    x1 = GradFS(KS*Hu(:,i));
    E(i) = F(x1) + G(K*x1);
end


sel = 1:round(options.niter/10);
clf;
loglog( E(sel)-min(E(:)) );
axis tight;