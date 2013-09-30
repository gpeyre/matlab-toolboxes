% test for impulsive noise removing
%
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');

if not(exist('name'))
    name = 'peppers';
    name = 'polygons_blurred';
    name = 'boat';
    name = 'barb';
    name = 'lena';
end

%% load image
n = 40;
n0 = [];
M = load_image(name, n0);
c = size(M)/2;
if strcmp(name, 'lena') && n==128
    c = [120 200]; % lena hat
end
if strcmp(name, 'mandrill') && n==128
    c = [350 350];
end
eta = 0; % .12;
M0 = rescale( crop(M,n,c), eta, 1-eta );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impulse noise
noise_type = 'impulse';
noise_type = 'laplacian';
switch noise_type
    case 'impulse'
        p = 0.05;
        M = compute_impulse_noise(M0,p);
    case 'laplacian'
        % exponent
        alpha = 1.6;
        % variance
        sigma = 0.07 * max( M0(:) );
        beta = 0;
        % deviation
        delta = 0;
        noise = stabrnd(alpha, beta, sigma, delta, n, n);
        noise = noise/std(noise(:)) * sigma;
        M = M0 + noise;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% median filter
k = 1;
disp('--> Median denoising.');
% Mmed = medfilt2(M, [w w], 'symmetric');
Mmed = perform_median_filtering(M,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NL median
disp('--> NLMedian denoising.');
options.max_dist = 6;
options.ndims = 25;
options.k = 3;
options.T = 0.04;
options.do_median = 1;

if 0
options.k = 1;
options.max_dist = k;
options.P = []; options.Psi = [];
options.Ha = ones(n+2*k);
options.H = ones(n+2*k);
Mnlmed = symmetric_extension(M,k);
end

Mnlmed = M;
[Mnlmed,Wx,Wy] = perform_nl_means(Mnlmed, options);

% Mnlmed = Mnlmed(k+1:end-k,k+1:end-k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NL mean
disp('--> NLMeans denoising.');
options.do_median = 0;
Mnlmean = M;
[Mnlmean,Wx,Wy] = perform_nl_means(Mnlmean, options);

pnoisy = psnr(M,M0);
pmed = psnr(M0,Mmed);
pnlmed= psnr(M0,Mnlmed);
pnlmean= psnr(M0,Mnlmean);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display
display_image_layout({M0 M Mmed Mnlmean Mnlmed}, ...
    {'Original' sprintf('Noisy,psnr=%.2f', pnoisy) sprintf('Median,psnr=%.2f', pmed) ...
        sprintf('NLmean,psnr=%.2f(+%.2f)', pnlmean, pnlmean-pmed) ...
        sprintf('NLmed,psnr=%.2f(+%.2f)', pnlmed, pnlmed-pmed) });