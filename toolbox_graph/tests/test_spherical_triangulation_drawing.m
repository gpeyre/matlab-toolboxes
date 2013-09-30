% test for drawing a spherical triangulation

rep = '../toolbox_graph_data/wrl/';
name = 'pawn';
name = 'test2';
name = 'test1';
[vertex,face] = read_wrl([rep name '.wrl']);

% project on sphere [should use here a spherical parameterization instead]
d = sqrt(sum(vertex.^2));
svertex = vertex./repmat(d, [3 1]);


% draw simply the projected mesh
clf;
plot_mesh(svertex,face);
shading faceted;

% now draw nice curve on the sphere
options.target_face = round(size(svertex,2)/2);
clf; 
plot_spherical_triangulation(svertex,face, options);