%   test for sub-images extraction from a target image
%   The set of images is obtained as small patches extracted from 
%   a big texture. The goal is to explore the structure of 
%   the local geometry of a texture.
%
%   Copyright (c) 2006 Gabriel Peyré

rep = '';
name = 'lena';
name = 'barb';

rep = 'data/textures/';
name = 'zigzag';
name = 'grass';
name = 'dunes';
name = 'thai_art';
name = 'warped_grid';
name = 'bricks2';
name = 'web';
name = 'turbulence1';
name = 'ecorce';
name = 'irreg_stones';
name = 'small_text';
name = 'sawtooth';
name = 'pasta';
name = 'sea_grid';
nbr = 1000;

w = 11; w1 = (w-1)/2;
options.dim = [1 1]*w;
options.sampling = 'rand';

M = load_image([rep name]);
M = sum(M,3);
options.equalize = 0;
X0 = load_subimages_dataset([rep name],nbr,options);


% turn it into a set of points
a = size(X0,1);
b = size(X0,2);
n = size(X0,3);
X = reshape(X0, a*b, n);

test_type = 'lle';
test_type = 'hlle';
test_type = 'geodesic';
test_type = 'isomap';
test_type = 'pca';

ndim = 3;

disp('-->Performing dimension reduction.');
switch lower(test_type)
    case 'isomap'
        clear options;
        options.nn_nbr = 8;
        xy = isomap(X, ndim, options);
    case 'lle'
        nn_nbr = 8;
        xy = lle(X,nn_nbr,ndim);
    case 'hlle'
        nn_nbr = 8;
        xy = hlle(X,nn_nbr,ndim);
    case 'pca'
        [tmp,xy] = pca(X,ndim);
    case 'geodesic'
        % take input points
        clf;
        imagesc(M); axis image; % axis off; 
        colormap gray(256);
        m = 6;
        [y,x] = ginput(m);  
        pos = round( [x';y'] );
        % add patches to list
        for i=1:size(pos,2)
            v = M(pos(1,i)-w1:pos(1,i)+w1,pos(2,i)-w1:pos(2,i)+w1);
            X(:,end+1) = v(:);
        end
        X0 = reshape(X,[a,b,size(X,2)]);
        options.landmarks = size(X,2)-m+1:size(X,2);
        % compute geodesic embeding
        xy = perform_geodesic_embedding(X,options);
end


if size(xy,1)>size(xy,2)
    xy = xy';
end


k = 30;
clf;
plot_flattened_dataset(xy,X0,k);
title(test_type);

return;

% perform exploration using an animation
clf;
hold on;
plot_scattered( xy(1:2,:) );
[y,x] = ginput( 2 );
hold off;
nbr_pts = 100;
path = [linspace(x(1),x(2),nbr_pts); linspace(y(1),y(2),nbr_pts)];
options.interp_mode = 'gaussian';
options.interp_mode = 'nearest';
A = perform_dimreduc_interpolation(path,xy(1:2,:),X,options);
A = reshape(A,[w w nbr_pts]);
play_movie(A);