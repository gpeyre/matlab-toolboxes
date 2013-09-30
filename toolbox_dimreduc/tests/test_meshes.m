% test for triangulation plotting
% (i.e. flattening of an disk-shaped 3D model)
%
%   Copyright (c) 2005 Gabriel Peyré

rep = 'data/';
name = 'nefertiti';
filename = [rep name '.off'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file loading
[vertex,face] = read_off(filename);
if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end
if size(face,1)<size(face,2)
    face = face';
end
A = triangulation2adjacency(face);
xy = vertex(:,1:2);

n = size(A,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original model
clf;
subplot(1,2,1);
plot_mesh(vertex,face);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use IsoMap

% build distance graph
D = build_euclidean_weight_matrix(A,vertex,Inf);
xy_isomap = isomap(D);

subplot(1,2,2);
gplot(A,xy_isomap,'k.-');
axis tight;
axis square;
axis off;
title('Isomap');