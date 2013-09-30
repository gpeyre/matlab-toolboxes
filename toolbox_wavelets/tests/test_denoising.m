% test for wavelet based denoising

name = 'lena';
rep = 'result_denoising/';

save_images = 1;
if save_images && exist(rep)~=7
    mkdir(rep);
end

n = 256;
M = load_image(name, n);
M = rescale(M, 10,245);

sigma = 15;
Mn = M + sigma * randn(n);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redundant wavelet transform

Jmin = 4;
options.wavelet_vm = 3;
options.wavelet_type = 'biorthogonal';
disp('Computing fwd transform.');
MW = perform_atrou_transform(Mn,Jmin,options);

% test several threshold and keep only best one
T_list = [2.5 3 3.5] * sigma;
err = [];
for T = T_list
    MWt = keep_above(MW, T);
    disp('Computing bwd transform.');
    M1 = perform_atrou_transform(MWt,Jmin,options);
    err = [err, psnr(M,M1)];
end

[tmp,I] = max(err);
T = T_list(I);
MWt = keep_above(MW, T);
disp('Computing bwd transform.');
M_wav_inv = perform_atrou_transform(MWt,Jmin,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orthogonal wavelet transform


Jmin = 4;
options.wavelet_vm = 4;
options.wavelet_type = 'biorthogonal';
disp('Computing fwd transform.');


use_rwt = 0;

if use_rwt
    h = daubcqf(4,'min');
    L = 3;
    [MW,L] = mdwt(Mn,h,L);
else
    MW = perform_wavelet_transform(Mn,Jmin,1,options);
end

% test several threshold and keep only best one
err = [];
for T = T_list
    MWt = keep_above(MW, T);
    disp('Computing bwd transform.');
    if use_rwt
        M1 = midwt(MWt,h,L);
    else
        M1 = perform_wavelet_transform(MWt,Jmin,-1,options);
    end
    err = [err, psnr(M,M1)];
end

[tmp,I] = max(err);
T = T_list(I);
MWt = keep_above(MW, T);
disp('Computing bwd transform.');
if use_rwt
    [M_wav_orth,L] = midwt(MWt,h,L);
else
    M_wav_orth = perform_wavelet_transform(MWt,Jmin,-1,options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering

% FFT-based wiener filtering (using the oracle fourier coefficients)
Mf = fft(Mn);
PMf = abs(Mf.*Mf); % power spectra
H = PMf./(PMf + n*sigma^2); % filter fourier transform
% compute convolution
Mnf = fft(Mn);
HMf = Mnf.*H;
M_wiener = real( ifft(HMf) );

% Matlab wiener filtering
disp('Computing Matlab Wiener filtering.');
M_wiener_matlab = wiener2(Mn,[5 5]);

% gaussian convolution 
err = [];
s = 4; % #pixels
h = compute_gaussian_filter([51 51],s/n, [n n]);
M_smooth = perform_convolution(Mn,h);



clf;
subplot(2,3,1);
imagesc(M);
axis image; axis off;
title('Original');

subplot(2,3,2);
imagesc(Mn);
axis image; axis off;
title( sprintf('Noisy (err=%.2fdB)', psnr(M,Mn)) );

subplot(2,3,3);
imagesc(M_wiener);
axis image; axis off;
title( sprintf('Wiener oracle Smoothing (err=%.2fdB)', psnr(M,M_wiener)) );

subplot(2,3,4);
imagesc(M_wiener_matlab);
axis image; axis off;
title( sprintf('Bloc-Wiener Denoised (err=%.2fdB)', psnr(M,M_wiener_matlab)) );


subplot(2,3,5);
imagesc(M_wav_orth);
axis image; axis off;
title( sprintf('Orthogonal Wavelet Denoised (T=%.2f, err=%.2fdB)', T, psnr(M,M_wav_orth)) );

subplot(2,3,6);
imagesc(M_wav_inv);
axis image; axis off;
title( sprintf('Redundant Wavelet Denoised (T=%.2f, err=%.2fdB)', T, psnr(M,M_wav_inv)) );

colormap gray(256);