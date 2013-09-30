% test for mesh simplification
name = 'cow.ascii.ply';
[vertex,face] = read_mesh(name);

nface = size(face,1);
rates = 0.1;

clf;
[vertex1,face1] = perform_mesh_simplification(vertex,face, round(nface*rates) );
plot_mesh(vertex,face);
figure;
plot_mesh(vertex1,face1);