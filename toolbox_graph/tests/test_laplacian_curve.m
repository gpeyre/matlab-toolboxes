%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for laplacian on a curve
% and try a comparison with the diffusion geometry
% basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = 'circle';
name = 'square';
n = 200;

% load curve
options.sampling = 'parabolic';
options.sampling = 'random';
options.sampling = 'uniform';
[c,t] = load_curve(name, n, options);

plot_curve(c, [], '.-');

% compute laplacian
clear options;
options.sigma = 0.1;
options.is_closed = 1;
options.method = 'combinatorial';
[Lc,Kc] = compute_laplacian_curve(c,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the eigenvalues, should be in [0,1], with decay 1/n
vc = sort(abs(eigs(Kc, n-1)));
clf;
plot( 1:n-1,vc(end:-1:1)  );
axis tight;
title('Eigenvalues of the diffusion kernel');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute basis function for laplacian
disp('Performing SVD.');
[Uc,S,V] = svd(full(Lc));
% reverse so that low frequencies comes first
Uc = Uc(:,end:-1:1);
% compute basis function for diffusion geometries
sigma = 0.1;
Options.Delta   = 1/sigma^2;
Options.NPoints = 2;
T = MakeDiffusion(c, 'Gauss', Options);
% T = Kc;
Ug = compute_diffusion_geometry(T, 1e-10);

% plot some eigenvectors
a = [2 5 10];
x = (1:n)';
clf;
subplot(2,1,1);
plot(x,Uc(:,a));
axis tight;
a = [2 5 10];
subplot(2,1,2);
plot(x,Ug(:,a));
axis tight;

% c = [(x>n/2) (x>n/2)];

% perform projection into laplacian basis
cc = Uc'*c;
cg = Ug'*c;
% error
ec = cumsum( cc(end:-1:1,1).^2 ) + cumsum( cc(end:-1:1,2).^2 );
ec = ec(end:-1:1);
eg = cumsum( cg(end:-1:1,1).^2 ) + cumsum( cg(end:-1:1,2).^2 );
eg = eg(end:-1:1);

% plot error
sel = 1:n/2;
clf;
loglog(sel,ec(sel), sel, eg(sel));
xlabel('log(#coefs)');
ylabel('log(error^2)');
axis tight;
legend('Laplacian','Diffusion geometries');
