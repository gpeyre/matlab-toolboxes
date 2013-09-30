% test for bilateral filtering

name = 'square';
name = 'lena';
n = 128;
M = load_image(name,n);

M = rescale(M,10,245);
sigma = 10;
Mn = M + sigma*randn(n);

% for limiting diffusion
sigma_list = [1 10 50];
% for scale of diffusion
winsize_list = [1 5 10];

clf;
for i=1:length(sigma_list)
    for j=1:length(winsize_list)
        sigma = sigma_list(i);
        winsize = winsize_list(j);
        M1 = perform_bilateral_filtering(Mn,winsize,sigma);
        subplot(length(sigma_list), length(winsize_list), i+(j-1)*length(sigma_list));
        imagesc( clamp(M1,0,255) ); 
        axis image; axis off;
        title( sprintf('spacial=%.1f  limiting=%.1f', winsize, sigma) );
    end
end
colormap gray(256);