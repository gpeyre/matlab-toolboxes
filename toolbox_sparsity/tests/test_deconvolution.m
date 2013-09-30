% test for sparsity-based deblurring

path(path, 'toolbox/');

n = 128;
name = 'lena';
sigma = 0.04; % noise level

% load the image
M0 = rescale( crop( load_image(name), n), 0.05 );

% options for the wavelet transform
Jmin = 4;
options.Jmin = Jmin;
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;

% width of the filter
options.eta = 3; 
% observation
Y = perform_blurring( M0, options.eta ) + sigma*randn(n);

% do inversion by iterative thresholding
lambda_min = 0;
lambda_max = 0.05;
nbr = 10;
lambda_list = linspace(lambda_max, lambda_min, nbr);

options.niter = 12; % iteration for inversion
options.mu = 1;
X = zeros(n^2,1); M = {}; perr = [];
for lambda = lambda_list
    options.lambda = lambda;
    options.X = X;
    X = perform_iterative_thresholding( @callback_blurring, Y(:), options );
    M{end+1} = perform_wavelet_transform(reshape(X,n,n), Jmin, -1, options);
    perr(end+1) = psnr(M0, M{end});
end
   
% plot error along regularization path
clf;
plot(lambda_list,perr);
xlabel('\lambda'); ylabel('PSNR');

% find best regularization lambda for inversion (this is cheating ...)
[pval,i] = max(perr);
clf;
imageplot(M0, 'Original', 1,3,1);
imageplot(Y, ['Observed, PSNR=' num2str(psnr(Y,M0))], 1,3,2);
imageplot(M{i}, ['Recovered, PSNR=' num2str(pval)], 1,3,3);

