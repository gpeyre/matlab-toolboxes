% test for block thresholding

name = 'peppers-bw';
name = 'barb';
name = 'lena';
name = 'boat';

n = 256;
c = [n n]/2;
if strcmp(name, 'lena') % && n==128
    c = [128, 220]; n = 128;
end

M0 = load_image(name);
M0 = rescale( crop(M0,n,c), .05,.95 );


%% add some noise
randn('state', 123456);
W = randn(n);
ptarget = 22; % boat
ptarget = 26; % boat
sigma = adjust_psnr(W,ptarget,M0);
M = M0 + sigma*W;

%% wavelets parameters
Jmax = log2(n)-1;
Jmin = 3;

options.wavelet_type = 'biorthogonal';
options.decomp_type = 'quad';
options.wavelet_vm = 4;
options.ti = 0;
options.block_size = 4;

if 0
MW = perform_wavelet_transform(M, Jmin, +1, options);
MWT = perform_thresholding(MW,T,thresh,options);
Mwav = perform_wavelet_transform(MWT, Jmin, -1, options);
end

type = 'waveortho';
if options.ti==1
    type = 'wavelet';
end

options.ti = 1;
options.thresh = 'hard';
[T,Mhard,pval,Tlist,err] = compute_best_threshold(type, M, M0, sigma, options);

options.ti = 0;
options.thresh = 'block';
[T,Mblock,pval,Tlist,err] = compute_best_threshold(type, M, M0, sigma, options);

pnoisy =    snr(M0,M);
pblock =      snr(M0,Mblock);
phard =      snr(M0,Mhard);

lgd = { ['Noisy, snr=' num2str(pnoisy)], ...
        ['Hard, snr=' num2str(phard)], ...
        ['Block, snr=' num2str(pblock)] };
imageplot({clamp(M) clamp(Mhard) clamp(Mblock)}, lgd, 2,2);

% disp(['SNR=' num2str(pwav)]);

