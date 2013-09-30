% tests for denoising using tv and wavelets

% 1D or 2D test
test = 1;

%% load data
if test==1
    name = 'piece-polynomial';
    n = 1024;
    x0 = load_signal(name, n);
    sigma = .03;
else
    name = 'barb';
    n = 256; 
    x0 = load_image(name);
    x0 = crop(x0,n);
    sigma = .05;
end
    
eta = .05;
x0 = rescale(x0,eta,1-eta);

randn('state',123456);
x = x0 + sigma*randn(size(x0));


%% wavelet domain denoising
%options for wavelets
options.ti = 1;
if test==1
    options.wavelet_vm = 4;
    options.wavelet_type = 'haar';
else
    options.wavelet_vm = 3;
    options.wavelet_type = 'biorthogonal';
    options.decomp_type = 'quad';
end
Jmin = 3;

niter = 10;
Tlist = sigma*linspace(.5,2.5,niter);
err = [];
if test==1
    xw = perform_wavelet_transform(x, Jmin, +1, options);
else
    xw = perform_atrou_transform(x, Jmin, options);
end
for i=1:niter
    xwT = perform_thresholding(xw,Tlist(i), 'soft');
    if test==1
        x1 = perform_wavelet_transform(xwT, Jmin, -1, options);
    else
        x1 = perform_atrou_transform(xwT, Jmin, options);
    end
    err(i) = norm(x0-x1,'fro');
    if i>1 && err(i)<min(err(1:i-1))
        xwav = x1;
    elseif i==1
        xwav = x1;
    end
end

%% TV denoising
options.niter = 5000;

options.etgt = norm(x-x0);
options.etgt = [];

if test==1
    options.lambda_max = .4;
    options.lambda_min = .1;
else
    options.lambda_max = .1;
    options.lambda_min = .0;
end
options.x0 = x0;
[xtv,err] = perform_tv_denoising(x,options);

pwav = psnr(x0,xwav,1);
ptv = psnr(x0,xtv,1);
pnoisy = psnr(x0,x1,1);
disp(['Wav=' num2str(pwav) 'dB, TV=' num2str(ptv) 'dB.' ]);

rep = 'results/tv-denoising/';
if not(exist(rep))
    mkdir(rep);
end

% save in a file the PSNR
filename = [rep 'results.txt'];
fid = fopen(filename, 'a');
fprintf( fid, '%s - noisy=%fdB - wav(%s)=%fdB - tv=%fdB\n', name, pnoisy, options.wavelet_type, pwav, ptv );
fclose(fid);


% display
if test==1
    fs = 20;
    clf;
    subplot(3,1,1);
    plot(x, 'k'); axis([1 n 0 1]);
    set(gca, 'FontSize', fs);
    subplot(3,1,2);
    plot(xwav, 'k'); axis([1 n 0 1]);
    set(gca, 'FontSize', fs);
    subplot(3,1,3);
    plot(xtv, 'k'); axis([1 n 0 1]);
    set(gca, 'FontSize', fs);
    saveas(gcf, [rep name '-tv-denoising.png'], 'png');
    saveas(gcf, [rep name '-tv-denoising.eps'], 'eps');
else
    imageplot({x0 x xwav xtv}, ...
            {'Original' ['Noisy=' num2str(pnoisy)] ['Wav=' num2str(pwav)]  ['TV=' num2str(ptv)] }, 2,2);
    saveas(gcf, [rep name '-tv-denoising.png'], 'png');
    warning off;
    imwrite(clamp(x0), [rep name '-original.png'], 'png' );
    imwrite(clamp(x), [rep name '-noisy.png'], 'png' );
    imwrite(clamp(xwav), [rep name '-wav.png'], 'png' );
    imwrite(clamp(xtv), [rep name '-tv.png'], 'png' );
    warning on;
end