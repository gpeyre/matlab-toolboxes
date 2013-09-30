% test for vector field reorientation

n = 128/2;
h = compute_gaussian_filter([61 61],20/(2*n),[n n]);
Phi = randn(n);
Phi = perform_convolution(Phi, h);
v0 = grad(Phi);


% corrupt
v = v0 .* repmat( sign(randn(n)), [1 1 2] );

% reorient
options.method = 'randomized';
options.niter_reorient = 500;
options.method = 'propagation';
options.method = 'laplacian';
v1 = v;
v1 = perform_vf_reorientation(v1, options);


clf;
imageplot(v0, 'original', 1,3,1);
imageplot(v, 'noisy', 1,3,2);
imageplot(v1, 'recovered', 1,3,3);