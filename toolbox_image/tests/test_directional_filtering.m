% test for directional smoothing

path(path, 'toolbox/');

n = 256;
name = 'lena';
M = load_image(name);
M = M(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2);


options.n_theta = 12;

% compute vector orthogonal to gradient
h = ones(7); h = h/sum(h(:));
Mh = perform_convolution(M,h);
V = compute_grad(Mh);
V = V(:,:,2:-1:1); V(:,:,1) = -V(:,:,1);

% perform directional filtering
options.sigma1 = 1;
options.sigma2 = 6;
M1 = perform_directional_filtering(M,V,options);

% perform isotropic filtering
options.sigma1 = 4;
options.sigma2 = 4;
M0 = perform_directional_filtering(M,V,options);

clf;
subplot(1,3,1);
imagesc(M); axis off; axis image;
title('Original');
subplot(1,3,2);
imagesc(M0); axis off; axis image;
title('Isotropic');
subplot(1,3,3);
imagesc(M1); axis off; axis image;
title('Anisotropic');
colormap gray(256);