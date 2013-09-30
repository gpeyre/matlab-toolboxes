% test for triangulation plotting
% (i.e. flattening of an disk-shaped 3D model)
%
%   Copyright (c) 2005 Gabriel Peyré

path(path, '../toolbox_graph_data/off/');

name = 'nefertiti';
name = 'mannequin';

rep = 'results/flattening/';
if not(exist(rep))
    mkdir(rep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file loading
[vertex,face] = read_off([name '.off']);

if strcmp(name, 'mannequin')
    vertex = -vertex;
end

if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end
if size(face,1)<size(face,2)
    face = face';
end
A = triangulation2adjacency(face);

n = length(A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original model
clf;
plot_mesh(vertex,face);
if strcmp(name, 'mannequin')
    view(28,10); camlight; 
end
% title('Original model');
saveas(gcf, [rep name '-mesh.png'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral graph drawing : use the eigenvectors of the laplacian
options.method = 'flattening';
options.laplacian = 'combinatorial'; 
xycombin = compute_parameterization(vertex,face, options);
clf;
plot_graph(A,xycombin);
% title('Combinatorial laplacian');
saveas(gcf, [rep name '-flattening-combinatorial.png'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same but use conformal laplacian
options.method = 'flattening';
options.laplacian = 'conformal'; 
xyconformal = compute_parameterization(vertex,face, options);
clf;
plot_graph(A,xyconformal);
% title('Conformal laplacian');
saveas(gcf, [rep name '-flattening-conformal.png'], 'png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use IsoMap
options.method = 'isomap';
xyisomap = compute_parameterization(vertex,face, options);
clf;
plot_graph(A,xyisomap);
% title('Isomap');
saveas(gcf, [rep name '-flattening-isomap.png'], 'png');


