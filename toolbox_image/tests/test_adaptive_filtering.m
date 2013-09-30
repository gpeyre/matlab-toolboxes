% test for adaptive filtering
%
% Perform a foveated filtering of the image.

path(path, 'toolbox/');

n = 256;
name = 'lena';
M = load_image(name, n);

% compute the filters
m = 31;
p = 30; % number of filters
sigma = linspace(0.05,6,p);
H = zeros(m,m,p);
for i=1:p
    H(:,:,i) = compute_gaussian_filter([m m],sigma(i)/n,[n n]);
end

% compute the index of the filters (for foveation)
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2 + Y.^2);
I = round(rescale(R,1,p));

tic;
M1 = perform_adaptive_filtering(M,H,I);
toc;

clf;
subplot(1,2,1);
imagesc(M); axis off; axis image;
subplot(1,2,2);
imagesc(M1); axis off; axis image;
colormap gray(256);