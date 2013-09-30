
path(path, '../toolbox_graph_data/off/');

name = 'mushroom';
name = 'nefertiti';
clear options;
options.name = name;
[vertex,face] = read_off([name '.off']);

rep = ['results/mesh-duality/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

clf;
plot_mesh(vertex,face,options);
shading faceted;
camlight; axis tight; axis equal;
saveas(gcf, [rep name '-mesh.png'], 'png');


clf;
A = triangulation2adjacency(face);
plot_graph(A,vertex);
view(2); axis tight; axis equal;
saveas(gcf, [rep name '-graph.png'], 'png');

[A1,vertex1] = compute_dual_graph(face,vertex);
clf;
plot_graph(A1,vertex1);
view(2); axis tight; axis equal;
saveas(gcf, [rep name '-dual.png'], 'png');

% 1:4 subdivision
options.sub_type = '1:4';
[vertex2,face2] = perform_mesh_subdivision(vertex',face',1, options);
[vertex3,face3] = perform_mesh_subdivision(vertex2,face2,1, options);

clf;
plot_mesh(vertex2,face2,options);
shading faceted;
camlight; axis tight; axis equal;
saveas(gcf, [rep name '-subdivide-4-once.png'], 'png');
clf;
plot_mesh(vertex3,face3,options);
shading faceted;
camlight; axis tight; axis equal;
saveas(gcf, [rep name '-subdivide-4-twice.png'], 'png');


% 1:3 subdivision
options.sub_type = '1:3';
[vertex2,face2] = perform_mesh_subdivision(vertex,face,1, options);
[vertex3,face3] = perform_mesh_subdivision(vertex2,face2,1, options);

clf;
plot_mesh(vertex2,face2,options);
shading faceted;
camlight; axis tight; axis equal;
saveas(gcf, [rep name '-subdivide-3-once.png'], 'png');
clf;
plot_mesh(vertex3,face3,options);
shading faceted;
camlight; axis tight; axis equal;
saveas(gcf, [rep name '-subdivide-3-twice.png'], 'png');
