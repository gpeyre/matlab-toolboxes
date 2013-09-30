% test for segmentation using a Linear/Non-linear/Linear (LNL) framework
% based on Gabor filtering
%
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');
path(path, 'images/');

if not(exist('M')) || not(exist('name'))
    names = {'fur', 'herringbone'};
    names = {'hair', 'reptilskin', 'fur', 'herringbone'};
    options.patchwork_mode = 'deterministic';
    options.patchwork_mode = 'random';
    %% load the image
    n0 = [];
    n = 220;
    name = [];
    M0 = {};
    nt = n/2; % size for training
    for i=1:length(names)
        M0{i} = load_image(names{i}, n0);
        M0{i} = rescale( sum(M0{i},3) );
        name = [name '-' names{i}];
    end
    name(1) = [];
    randn('seed',666);
    [M,B0] = compute_texture_patchwork(M0,n,options);
end
ntextures = length(M0);


rep = 'results/segmentation-gabor/';
if not(exist(rep))
    mkdir(rep);
end


%% compute a set of gabor features
options.gabor_mode = 'radial';
options.gabor_mode = 'oriented';
options.iscomplex = 1;
options.ntheta = 8;
options.nsigma = 6;
options.nfreq = 3;
options.add_spacial = 0;
fprintf('--> Computing Gabor filterings: ');
[E,F] = compute_gabor_features(M,options);


%% perform various blurring and find best segmentation
fprintf('--> Test various blurring width: ');
options.ntextures = ntextures;
options.segmentation_method = 'kmeans';
options.oracle = B0;
[B,err] = perform_segmentation(E,options);

%% colorize the segmentation
Bcol = perform_segmentation_colorization(B);
B0col = perform_segmentation_colorization(B0);

%% display
clf;
imageplot({cat(3,M,M,M),B0col,Bcol}, {'Input', 'Ground trust', ['Segmentation ' num2str(err*100,3) '%']});



clf;
subplot(1,2,1);
display_segmentation(B0,M);
title('Ground trust');
subplot(1,2,2);
display_segmentation(B,M);
title(['Segmentation ' num2str(err*100,3) '%']);


saveas(gcf, [rep name '-segmentation-gabor.png'], 'png');
% nt' num2str(ntheta) '-ns' num2str(nsigma) '-spa'  num2str(add_spacial)
% '-' gabor_mode '