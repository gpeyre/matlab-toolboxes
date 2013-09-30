% test for tetrahedral mesh loading and display
path(path,'../toolbox_graph_data/tet/');

name = 'hand';
name = 'skull';
name = 'cow';
name = 'arma00';
name = 'torso';

[vertex,face] = read_tet([name '.tet']);
n = size(vertex,2);

options.cutting_interactive = 1;
plot_mesh(vertex,face);