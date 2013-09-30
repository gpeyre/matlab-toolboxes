% test for gaussian filtering
n = 100;
M = [ones(n/2,n);zeros(n/2,n)];
M = M + 0.1*randn(n);   % adding some noise

subplot(2,2,1);
imagesc(M);
title('Original Image');
axis image; axis off;

% without convolution
[gx,gy] = gradient(M);
subplot(2,2,2);
imagesc(gx.^2 + gy.^2);
title('Gradient');
axis image; axis off;

% with convolution
h = compute_gaussian_filter([20 20],0.02,n);
subplot(2,2,3);
imagesc(h);
title('Filter');
axis image; axis off;

Mh = perform_convolution(M,h);
[gx,gy] = gradient(Mh);
subplot(2,2,4);
imagesc(gx.^2 + gy.^2);
title('Smoothed gradient');
axis image; axis off;

colormap gray(256);