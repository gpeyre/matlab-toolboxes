% test for denoising



name = 'lena';
name = 'peppers';
name = 'boat';
name = 'lena';
name = 'polygons_blurred';
name = 'barb';

n = 128*2;
M = load_image(name);
eta = 0.12;
c = size(M)/2;
% c = [120 200]; % lena hat
M0 = rescale( M(c(1)-n/2+1:c(1)+n/2,c(2)-n/2+1:c(2)+n/2 ), eta, 1-eta );

% add some noise
sigma = 0.04 * max( M0(:) );
randn('state',1234);    % to have reproductible results
M = M0 + randn(n)*sigma;

options.sigma = sigma;
options.repres2 = 'daub3';
options.parent = 1;
M1 = perform_blsgsm_denoising(M, options);

pnoisy = psnr(M,M0);
pwav = psnr(M0,M1);

display_image_layout({M0 M M1}, ...
    {'Original' sprintf('Noisy,psnr=%.2f', pnoisy) sprintf('Denoised,psnr=%.2f', pwav)  });