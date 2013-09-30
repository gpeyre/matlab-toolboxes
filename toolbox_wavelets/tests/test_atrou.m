% test for a trou transform
n = 256;
name = 'lena';

M = load_image(name, n);

Jmin = 3;
options.wavelet_vm = 3;
options.wavelet_type = 'biorthogonal';
options.decomp_type = 'quad';
MS = perform_atrou_transform(M,Jmin,options);

rep = ['results/' name '/'];

save_images = 1;
if exist(rep)~=7
    mkdir(rep);
end

warning off;
for s=1:9
    imwrite(rescale(MS{s}), [rep name '_wavtrou_' num2str(s) '.png'], 'png');
end
warning on;