% test NL means denoising
%
%   Copyright (c) 2007 Gabriel Peyre

path(path,'toolbox/');
path(path, 'images/');

%% load the image
name = 'lenacoul';
name = 'barb';
m = 50; n = 90;     % just to show that it works with rectangular images
n = m;
M = load_image(name);
% crop the image
M = rescale( crop(M, [m n]) );

sigma = 0.03; % variance of additional noise
if sigma>0
    % avoid saturation
    M = clamp( rescale(M,sigma,1-sigma) + sigma * randn(size(M)) );
end


%% options of NL means
options.k = 3;          % half size for the windows
options.T = 0.03;       % width of the gaussian, relative to max(M(:))  (=1 here)
options.max_dist = 15;  % search width, the smaller the faster the algorithm will be
options.ndims = 30;     % number of dimension used for distance computation (PCA dim.reduc. to speed up)
options.do_patchwise = 0;

%% do denoising
tic;
[M1,Wx,Wy] = perform_nl_means(M, options);
toc;

clf;
imagesc(M1);
colormap gray(256);
return;

%% display results
ax = [];
clf;
ax(1) = subplot(2,2,1);
imagesc(M); axis image; axis off;
title('Original image');
ax(2) = subplot(2,2,2);
imagesc(clamp(M1)); axis image; axis off;
title('Denoised');
ax(3) = subplot(2,2,3);
imagesc(rescale(M-M1));
title('Removed noise');
axis image; axis off;
if size(M,3)>1
    colormap gray(256);
end
linkaxes(ax,'xy');