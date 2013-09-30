% test for the resolution of anlysis-based regularization using Chambolle
% algorithm. The goal is to solve, given f of length n, for
%       min_g  1/2*|f-g|^2 + lambda * |A*g|_1
% where A is an (m,n) operator and |x|_1=\sum_k |x(k)|

path(path, 'toolbox/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first a 1D example using derivatives
n = 512;
f0 = rescale( load_signal('Piece-Regular', n) );
sigma = 0.04;
f = f0+sigma*randn(n,1);

e = ones(n,1);
G = spdiags([e -e], [0 1], n, n);
G(end,:) = [];

options.niter = 150;
options.eta = 0.1;
options.lambda = 0.3;
options.lambda_min = 0.001; 
options.lambda_max = 0.15;
[g,g_list,E] = perform_analysis_regularization(f, G, options);

% select the optimal threshold value by oracle
err = sum( (g_list - repmat(f0,[1 options.niter])).^2 );
[tmp,k] = min(err);

clf;
subplot(3,1,1);
plot(f0); axis tight; title('Original');
subplot(3,1,2);
plot(f); axis tight; title('Noisy');
subplot(3,1,3);
plot(g_list(:,k)); axis tight; title('Regularized');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D example using TI wavelets
name = 'lena';
n = 256;
M0 = rescale( crop( load_image(name), n ), 0.05,0.95 );
sigma = 0.1;
M = M0+sigma*randn(n);


options.niter = 40;
options.eta = 0.1;
options.lambda_min = 0.01; 
options.lambda_max = 0.13;
options.n = n;

[G,G_list] = perform_analysis_regularization(M(:), @callback_ti_wavelets, options);

% select the optimal threshold value by oracle
err = sum( (G_list - repmat(M0(:),[1 options.niter])).^2 );
[tmp,k] = min(err);
G = reshape( G_list(:,k), n,n );

% simple denoising by thresholding at 3sigma
MW = callback_ti_wavelets(M(:),+1,options);
sigma_list = sigma*linspace(2.5,3,12);
nsigma = length(sigma_list);
H_list = zeros(n,n,nsigma);
for i=1:nsigma
    progressbar(i,nsigma);
    MWT = perform_thresholding(MW,sigma_list(i),'hard');
    H_list(:,:,i) = reshape( callback_ti_wavelets(MWT,-1,options), n,n);
end
err = sum( sum( (H_list - repmat(M0,[1 1 nsigma]) ).^2, 1 ), 2);
[tmp,k] = min(err);
H = H_list(:,:,k);

clf;
imageplot(M0, 'Original', 2,2,1);
imageplot(clamp(M), ['Noisy, psnr=' num2str(psnr(M,M0),3)], 2,2,2);
imageplot(clamp(G), ['Synth.denoising, psnr=' num2str(psnr(G,M0),3)], 2,2,3);
imageplot(clamp(H), ['Thresh.denoising, psnr=' num2str(psnr(H,M0),3)], 2,2,4);
