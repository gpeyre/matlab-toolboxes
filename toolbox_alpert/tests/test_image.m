% test for image transformation

name = 'barb';
n = 128;

M = load_image(name, n);

% find sampling location for an image
[Y,X] = meshgrid(1:n,1:n); 
pos = [X(:),Y(:)]';


% perform the transform
alpert_vm = [2 2];
% this is very important for images
options.part_type = '2axis';  
w = perform_alpert_transform_2d(M(:),pos,alpert_vm, 1, options);

% perform some thresholding
T = 20; % this threshold is ok for images in range [1...256], otherwise scale it
w1 = w .* (abs(w)>T);

% reconstruction
v1 = perform_alpert_transform_2d(w1,pos,alpert_vm, -1, options);
M1 = reshape(v1,n,n);

clf;
subplot(1,2,1);
imagesc(M);
title('Original'); axis image; axis off;
subplot(1,2,2);
imagesc(M1);
title('Approximated'); axis image; axis off;
colormap gray(256);