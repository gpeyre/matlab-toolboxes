% test for denoising using Non Local Means vs Wavelets
%
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');

if not(exist('name'))
    name = 'peppers';
    name = 'polygons_blurred';
    name = 'boat';
    name = 'mandrill';
    name = 'lenacoul';
    name = 'lena';
    name = 'images/corral';
    name = 'barb';
end

%% load image
n = 64*2;
n0 = [];
M = load_image(name, n0);
c = size(M)/2;
if strcmp(name, 'lena') && n==128
    c = [120 200]; % lena hat
end
if strcmp(name, 'mandrill') && n==128
    c = [350 350];
end
eta = 0.12;
M0 = rescale( crop(M,n,c), eta, 1-eta );
s = size(M0,3); % number of colors

%% add some Gaussian noise
sigma = 0.06 * max( M0(:) );
randn('state',1234);    % to have reproductible results
M = M0 + randn(n,n,s)*sigma;
% approx wavelet threshold
T = 3*sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TI Wavelets
clear options;
options.wavelet_type = 'biorthogonal';
options.wavelet_vm = 3;
Jmin = 4;
options.decomp_type = 'quad';
MW = perform_atrou_transform(M,Jmin,options);
disp('--> Computing best threshold (wavelets).');
Twav = compute_best_threshold('wavelet',M,M0,sigma, options);
MWT = perform_thresholding(MW, T);
MW1 = perform_atrou_transform(MWT,Jmin,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Portilla and Simoncelli method
options.sigma = sigma;
options.repres2 = 'daub3';
options.parent = 1;
disp('--> BLS-GSM denoising.');
ML1 = perform_blsgsm_denoising(M, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-local means
% options of NL means  (15,30)
options.max_dist = 5;  % search width, the smaller the faster the algorithm will be
options.ndims = 40;     % number of dimension used for distance computation (PCA dim.reduc. to speed up)
options.do_patchwise = 0;
options.mask = 'linear';
options.mask = 'cst';
disp('--> NLMeans denoising.');
options.Tlist = linspace(0.02,0.07,12);
options.klist = [4 5];
[tmp,options] = compute_best_threshold('nlmeans',M,M0,sigma, options);
[MN1,Wx,Wy] = perform_nl_means(M, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-local patchwise
do_patchwise = 1;
if do_patchwise
options.do_patchwise = 1;
disp('--> NLMeans denoising.');
[tmp,options] = compute_best_threshold('nlmeans',M,M0,sigma, options);
[MNpwise1,Wx,Wy] = perform_nl_means(M, options);
end



pnoisy = psnr(M,M0);
pwav = psnr(M0,MW1);
pbls= psnr(M0,ML1);
pnlm= psnr(M0,MN1);
if do_patchwise
    pnlpwise = psnr(M0,MNpwise1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display
imgs = {M0 M MW1 ML1 MN1};
opts =  {'Original' sprintf('Noisy,psnr=%.2f', pnoisy) sprintf('Wav,psnr=%.2f', pwav) ...
        sprintf('BLS,psnr=%.2f(+%.2f)', pbls, pbls-pwav) ...
        sprintf('NLMeans,psnr=%.2f(+%.2f)', pnlm, pnlm-pwav) };
if do_patchwise
    imgs{end+1} = MNpwise1;
    opts{end+1} = sprintf('NLMedian,psnr=%.2f(+%.2f)', pnlpwise, pnlpwise-pwav);
end
display_image_layout(imgs, opts );

% save
repimg = ['results/denoising/'];
if not(exist(repimg))
    mkdir(repimg);
end
saveas(gcf, [repimg name '-denoising.png'], 'png');

repimg = ['results/denoising/' name '/'];
if not(exist(repimg))
    mkdir(repimg);
end


%% write results
fid = fopen([repimg name '-results.txt'], 'a');
fprintf(fid, '--> %s: n=%d, max_dist=%d, ndims=%d, sigma=%.2f, k=%d, T=%.2f, mask=%s\n', name, n, ...
        options.max_dist, options.ndims, sigma, options.k, options.T, options.mask);
fprintf(fid, 'Noisy    %.4f\n', pnoisy);
fprintf(fid, 'Wavelets %.4f\n', pwav);
fprintf(fid, 'BLS-GSM  %.4f (+%.4f)\n', pbls, pbls-pwav);
fprintf(fid, 'NL-Means %.4f (+%.4f)\n', pnlm, pnlm-pwav);
if do_patchwise
    fprintf(fid, 'NL-Median %.4f (+%.4f)\n', pnlpwise, pnlpwise-pwav);
end
fclose(fid);

% scaling for method noise
s = std(M0(:)-MW1(:))*6;

%% save images
warning off;
imwrite(clamp(M0), [repimg name '-original.png'], 'png');
imwrite(clamp(M), [repimg name '-noisy.png'], 'png');
imwrite(clamp(MW1), [repimg name '-wavelets.png'], 'png');
imwrite(clamp(ML1), [repimg name '-blsgsm.png'], 'png');
imwrite(clamp(MN1), [repimg name '-nlmeans.png'], 'png');
imwrite(clamp( 0.5+ (M-MW1)/s ) , [repimg name '-wavelets-methnoise.png'], 'png');
imwrite(clamp( 0.5+ (M-ML1)/s ) , [repimg name '-blsgsm-methnoise.png'], 'png');
imwrite(clamp( 0.5+ (M-MN1)/s ) , [repimg name '-nlmeans-methnoise.png'], 'png');
if do_patchwise
imwrite(clamp(MNpwise1), [repimg name '-nlpatchwise.png'], 'png');    
imwrite(clamp( 0.5+ (M-MNpwise1)/s ) , [repimg name '-nlpatchwise-methnoise.png'], 'png');
end
warning off;