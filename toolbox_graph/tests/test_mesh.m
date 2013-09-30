% test for the various functionality for mesh processing

%% load the mesh
name = 'mushroom';
[vertex,face] = read_mesh([name '.off']);

%% display
clf;
plot_mesh(vertex, face);

%% compute the adjacency matrix
A = triangulation2adjacency(face);

%% display the mesh as a graph
clf;
plot_graph(A,vertex);

%% compute the dual graph
[A1,vertex1] = compute_dual_graph(face,vertex);
clf;
plot_graph(A1,vertex1);

%% compute a graph that links center of faces to vertex
[Ad,vertexc] = compute_semidual_graph(face, vertex);
clf;
plot_graph(Ad,vertexc);

%% vertex ring
vring = compute_vertex_ring(face);
%% face ring
fring = compute_face_ring(face);
