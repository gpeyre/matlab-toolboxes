% test for isomap with landmarks points

% simple tests for isomap/lle on 3D data
%
%   Copyright (c) 2005 Gabriel Peyré


name = 'square';
name = 'scurve';
name = '3d_cluster';
name = 'puncted_sphere';
name = 'swissroll';


n = 1000;
options.sampling = 'uniform';
options.sampling = 'rand';
options.noise_level = 0;
[X,col] = load_points_set( name, n, options );

% number of neighbors
nn_nbr = 6;

clf;
subplot(1,3,1);
plot_scattered(X,col);
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% without landmarks
clear options;
options.nn_nbr = nn_nbr;
xy_isomap = isomap(X,2, options);

subplot(1,3,2);
plot_scattered(xy_isomap,col);
axis off;
title('Isomap');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with landmarks
p = 20; % number of ldm
landmarks = unique( floor(rand(p,1)*n)+1 );

options.landmarks = landmarks;
xy_landmarks = isomap(X,2, options);

subplot(1,3,3);
hold on;
plot_scattered(xy_landmarks,col);
plot(xy_landmarks(landmarks,1),xy_landmarks(landmarks,2), '*');
hold off;
axis off;
title('Isomap with Landmarks');
