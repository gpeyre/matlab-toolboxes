% test for vector quantization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '3d';
name = 'lena';

if strcmp( name, 'lena' )
    n = 128;
    X = load_image(name, n);
    A = [2,2];
else
    n = 10;
    x = linspace(0,1,n);
    [X,Y,Z] = meshgrid(x,x,x);
    % some smooth signal
    X = X.^2+Y.^2+Z.^2;
    A = [2 1 1];
end



options.block_size = A;
options.nbr_iterations = 20;

% number of quantized vector
p = prod(size(X)) / prod(A);

options.nbr_vectors = round(0.2*p);

% vector used to construct dictionnary
options.nbr_samples = max( round(p/2), options.nbr_vectors/10 );
options.nbr_samples = 1000;

[C,IDX] = perform_vector_quantization(X,options);
X1 = perform_vector_quantization(C,IDX,options);

% reconstruction error
err = psnr(X,X1);
% number of bits
b = log2( size(C,2) ); % number of bits per index
N = prod(size(X));
bpp_idx = prod(size(IDX))*b / N;
bpp_c = prod(size(C))*8 / N;

% plot the two images
clf;
subplot(1,2,1);
imagesc(X(:,:,1));
axis image; axis off; colormap gray(256);
subplot(1,2,2);
imagesc(X1(:,:,1));
axis image; axis off; colormap gray(256);
title( sprintf('IDX=%.2fbpp C=%.2fbpp TOT=%.2f  psnr=%.2fdB', bpp_idx, bpp_c, bpp_idx + bpp_c, err) );