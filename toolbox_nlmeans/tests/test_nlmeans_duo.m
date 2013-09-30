% test NL means denoising using two different images
%
%   Copyright (c) 2007 Gabriel Peyre


%% load the image
rep = 'images/';    % directory where you can find the image
rep = '../../images/';

name = 'lenacoul';
name = 'monet';
name = 'hair';
name = 'corral';
name = 'dead_leaf';
name = 'reptilskin';
name = 'warped_grid';
name = 'simoncelli7';

use_same_image = 1;

m = 60; n = 80;
m = 128; n = m;
M = load_image([rep name]);
% crop the image
M = rescale( M(end/2-m/2+1:end/2+m/2,end/2-n/2+1:end/2+n/2,:) );

%% you can use another image for denoising
name1 = 'olives';
name1 = 'hair';
ma = 128; na = 65; % to test different size
ma = m; na = n;
Ma = load_image([rep name1]);
Ma = rescale( Ma(end/2-ma/2+1:end/2+ma/2,end/2-na/2+1:end/2+na/2,:) );

if use_same_image
    M00 = load_image([rep name]);
    M00 = rescale( sum(M00,3) );
    M0 = rescale(M00,0.05,0.95);
    clf;
    imagesc(M0); axis image; axis off;
    colormap gray(256);
    [y,x] = ginput(2);
    x = clamp(x,n/2,size(M0,1)-n/2);
    y = clamp(y,n/2,size(M0,1)-n/2);
    Ma = rescale( crop(M0, na, [x(2), y(2)]) );
    M0 = rescale( crop(M0, n, [x(1), y(1)]) );
    
    % Ma = perform_histogram_equalization(Ma,M0);
else
    %% to enhanced the results, we first match the colors
    Ma = perform_histogram_equalization(Ma,M);
end

clf;
subplot(1,2,1);
imagesc(M0); axis image; axis off;
title('Orignal');
subplot(1,2,2);
imagesc(Ma); axis image; axis off;
title('Target');
if size(M0,3)==1
    colormap gray(256);
end

%% options of NL means
options.k = 2;          % half size for the windows
options.max_dist = 30;  % search width
options.ndims = 25; % number of dimension used for distance computation (PCA dim.reduc.)

%% add some noise
sigma = 0.08*max(abs(M(:)));
M = M0 + randn(size(M0))*sigma;

%% test for a wide range of kernel widths
err_duo = [];
err_mono = [];
niter = 8;
Tlist = linspace(0.5, 1.2, niter) * sigma;

%% do denoising, my recommandation is to do several iterations
for i=1:niter
    progressbar(i,niter);
    T = Tlist(i);
    options.T = T;       % width of the gaussian, relative to max(M(:))  (=1 here)
    options.Ma =  [];
    M_mono = perform_nl_means(M, options);
    options.Ma =  Ma;
    M_duo = perform_nl_means(M, options);
    err_mono(end+1) = psnr(M0,M_mono);
    err_duo(end+1) = psnr(M0,M_duo);
end

clf;
plot(Tlist/sigma,err_mono, Tlist/sigma,err_duo);
legend('mono', 'duo'); axis tight;
saveas(gcf, [rep name '-psnr.png'], 'png');

%% find min-error
[pmono,i] = max(err_mono);
Tmono = Tlist(i);
[pduo,i] = max(err_duo);
Tduo = Tlist(i);
pnoisy = psnr(M,M0);

%% do final denoising
options.T = Tmono;       % width of the gaussian, relative to max(M(:))  (=1 here)
options.Ma =  [];
M_mono = perform_nl_means(M, options);
options.T = Tduo; 
options.Ma =  Ma;
M_duo = perform_nl_means(M, options);

rep = ['results/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

%% display results
clf;
subplot(2,2,1);
imagesc(clamp(Ma)); axis image; axis off;
title('Source');

subplot(2,2,2);
imagesc(clamp(M)); axis image; axis off;
title(['Noisy, PSNR=' num2str(pnoisy,4)]);

subplot(2,2,3);
imagesc(clamp(M_mono)); axis image; axis off;
title(['Denoised mono, PSNR=' num2str(pmono,4)]);

subplot(2,2,4);
imagesc(clamp(M_duo)); axis image; axis off;
title(['Denoised duo, PSNR=' num2str(pduo,4)]);

if size(M,3)==1
    colormap gray(256);
end
saveas(gcf, [rep name '-result.png'], 'png');

warning off;
imwrite(clamp(M00), [name '.png'], 'png');
imwrite(clamp(Ma), [name '-src.png'], 'png');
imwrite(clamp(M0), [name '-orig.png'], 'png');
imwrite(clamp(M), [name '-noisy.png'], 'png');
imwrite(clamp(M_mono), [name '-mono.png'], 'png');
imwrite(clamp(M_duo), [name '-duo.png'], 'png');
warning off;



