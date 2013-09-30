% fill missing part of 3D models

rep = '../toolbox_graph_data/wrl/';
path(path, rep);

name = 'venus';
[vertex,face] = read_mesh(name);

boundary=compute_boundary(face);

clf;
hold on;
plot_mesh(vertex,face);
shading faceted;
h = plot3(vertex(1,boundary),vertex(2,boundary),vertex(3,boundary), 'r.');
set(h,'MarkerSize', 20);

% perform pca projection
% [tmp,pts] = pca(vertex(:,boundary),2);

vertex(:,end+1) = mean(vertex(:,boundary)')';
boundary(end+1) = size(vertex,2);

m = length(boundary)-1;
t = linspace(0,2*pi,m+1); t(end) = [];
pts = [cos(t); sin(t)];
pts(:,end+1) = [0;0];

% triangulate
tri = delaunay(pts(1,:),pts(2,:))';

% check wether the triangles are in the correct orientation


clf;
hold on;
plot_mesh(pts,tri);
shading faceted;
h = plot(pts(1,:),pts(2,:), '.r');
set(h,'MarkerSize', 20);

face = [face, boundary(tri)];

clf;
plot_mesh(vertex,face);
shading faceted;

write_wrl([rep name '-filled.wrl'], vertex,face);