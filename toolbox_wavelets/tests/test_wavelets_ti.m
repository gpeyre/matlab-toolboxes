%% Test for several implementation of wavelet transform

Jmin = 3;
options.ti = 1;
options.wavelet_vm = 3;
for dimension=1:2

    if dimension==1
        x = load_signal('Piece-Regular', 512);
    else
        x = load_image('lena', 256);
    end

    %% test for RWT
    options.use_mex = 1;
    options.wavelet_type = 'daubechies';
    tic;
    y = perform_wavelet_transform(x, Jmin, +1, options);
    x1 = perform_wavelet_transform(y, Jmin, -1, options);
    disp(['RWT,     Time=' num2str(toc) ', Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);

    %% test for wavelab
    options.use_mex = 0;
    options.wavelet_type = 'daubechies';
    tic;
    y = perform_wavelet_transform(x, Jmin, +1, options);
    x1 = perform_wavelet_transform(y, Jmin, -1, options);
    disp(['Wavelab, Time=' num2str(toc) ', Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);
    
    
    %% test for lifting
    options.use_mex = 0;
    options.wavelet_type = 'biorthogonal';
    tic;
    y = perform_wavelet_transform(x, Jmin, +1, options);
    x1 = perform_wavelet_transform(y, Jmin, -1, options);
    disp(['Lifting, Time=' num2str(toc) ', Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);

    %% test for LIW
    if dimension==2
        options.use_mex = 1;
        options.wavelet_type = 'biorthogonal';
        tic;
        y = perform_wavelet_transform(x, Jmin, +1, options);
        x1 = perform_wavelet_transform(y, Jmin, -1, options);
        disp(['LIW,     Time=' num2str(toc) ', Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);
    end
end

%% compare LIW vs. Lifting for denoising XP
M0 = load_image('lena');
M0 = crop(M0,256);
M0 = rescale(M0,.03,.97);
sigma = .1;
M = M0 + randn(size(M0))*sigma;

Jmin = 3;
options.ti = 1;
options.wavelet_vm = 3;
options.wavelet_type = 'biorthogonal';

for i = 0:1
    options.use_mex = i;
    tic;
    MW{i+1} = perform_wavelet_transform(M,Jmin,+1,options);
    MW{i+1} = perform_thresholding(MW{i+1},3*sigma);
    Mres{i+1} = perform_wavelet_transform(MW{i+1},Jmin,-1,options);
    toc;
end
lgd = { num2str(snr(M0,Mres{1})) num2str(snr(M0,Mres{2})) };
imageplot(Mres,lgd);