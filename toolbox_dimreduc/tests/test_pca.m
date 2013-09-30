% test for pca versus empca

p = 1000; % # points
d = 200; % number of dimension
X = rand(d,p);
numvecs = 20;

% classical eigs decompositon
tic;
options.use_em = 0;
[Y,X1,v,Psi] = pca(X,numvecs, options);
toc;

% EM algorithm
tic;
options.use_em = 1;
options.iter = 20;
[Y,X1,v,Psi] = pca(X,numvecs, options);
toc;