% Demo of the Laplacian pyramid functions
n = 256;
M = load_image('lena', n);

% Laplacian decomposition using 9/7 filters and 5 levels
pfilt = '9/7'; J = 5;
MP = perform_pyramid_transform_do(M, pfilt, J);
% Wavelet transform
Jmin = 4;
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;
MW = perform_wavelet_transform(M, Jmin, +1, options);


% Display output of the pyramid
figure(1)
colormap(gray);
nr = floor(sqrt(J+1));
nc = ceil((J+1)/nr);
for l = 1:J+1
    subplot(nr, nc, l); 
    imagesc(MP{l});
    axis image; axis off;
end

% Reconstruction
m = floor( 0.05*n^2 );
MPt = keep_biggest(MP, m);
MWt = keep_biggest(MW, m);
M1 = perform_pyramid_transform_do(MPt, pfilt);
M2 = perform_wavelet_transform(MWt, Jmin, -1, options);

% Show perfect reconstruction
figure(2);
colormap gray;
subplot(1,3,1), imagesc(M);
axis image; axis off;
title('Original image');
subplot(1,3,2), imagesc(M1);
axis image; axis off;
title(sprintf('Pyramid, PSNR=%.2f dB', PSNR(M, M1)))
subplot(1,3,3), imagesc(M2);
axis image; axis off;
title(sprintf('Wavelet, PSNR=%.2f dB', PSNR(M, M2)))