% test for line integral convolution

n = 180;
path(path, 'toolbox/');

rep = 'results/lic/';
if not(exist(rep))
    mkdir(rep);
end

options.bound = 'per';
%% generate a random irrotational vector field
sigma = 30;
v = perform_blurring(randn(n,n,2), sigma, options);
v = perform_vf_normalization(v);

options.histogram = 'gaussian';
options.histogram = 'linear';
options.verb = 1;
% size of the features
options.spot_size = 1.3;


%% original image
name = 'mandrill-color';
name = 'flowers';
name = 'rand';
options.dt = 0.5;
w_list = [4 6 8 10 12 14];
if strcmp(name, 'rand')
    sigma = 1.2;
    M0 = perform_blurring(randn(n), sigma, options);
    M0 = perform_histogram_equalization(M0, options.histogram);
else
    options.dt = 1.5;
    w_list = w_list * 2;
    options.histogram = [];
    M0 = load_image(name,n);
    M0 = rescale( crop(M0,n) );
end

%% iterated lic
M = M0;
w = 12;
for i=1:4
    options.M0 = M;
    M = perform_lic(v, w, options);
end


%% lic for increasing times
close all;
clf;
options.M0 = M0;
for i=1:min(length(w_list),6)
    M = perform_lic(v, w_list(i), options);
    imageplot(M, '', 2,3, i);
end

saveas(gcf, [rep name '-lic-results.png'], 'png');


return;

%% display with overlapping vector field
if size(M,3)==1
    figure; clf;
    sub = 4;
    plot_vf(v(1:sub:end, 1:sub:end,:), M);
    colormap gray(256);
    saveas(gcf, [rep, 'lic-vector-field.png'], 'png');
end