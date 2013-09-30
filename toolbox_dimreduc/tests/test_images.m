% test for dimension reduction on images datasets

name = 'binaryalphadigs';
name = 'umist_cropped';
name = 'edges';
name = 'disks';
name = 'frey_rawface';
name = 'digits';
name = 'olivettifaces';
name = 'teapots-bw';

options.nclass = 4;
options.nbr = 500;

if strcmp(name, 'umist_cropped')
    options.nclass = 1:5;
end
if strcmp(name, 'disks')
    options.smoothing = 0.05;
    options.dim = [64 64];
    options.nbr = 500;
end
if strcmp(name, 'edges')
    options.smoothing = 0;
    options.dim = [20 20];
    options.nbr = 25;
end

disp('-->Reading database');
M = load_images_dataset(name, options);

% turn it into a set of points
a = size(M,1);
b = size(M,2);
n = size(M,3);
X = reshape(M, a*b, n);

nn_nbr = 5;

% test_type = 'hlle';
test_type = 'lle';
test_type = 'pca';
test_type = 'isomap';

disp('-->Performing dimension reduction.');
if strcmp(test_type, 'isomap')
    clear options;
    options.nn_nbr = nn_nbr;
    xy = isomap(X,2, options);
elseif strcmp(test_type, 'lle')
    nn_nbr = 16;
    xy = lle(X,nn_nbr,2);
elseif strcmp(test_type, 'hlle')
    xy = hlle(X,nn_nbr,2);
elseif strcmp(test_type, 'pca')
    [tmp,xy] = pca(X,2); xy = xy';
else
    error('Unknown method.');
end


k = 30;
clf;
plot_flattened_dataset(xy,M,k);
title(test_type);