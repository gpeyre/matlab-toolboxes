%   test for triangulation parameterization
%   using free boundariy (neumann conditions)
%
%  Copyright (c) 2007 Gabriel Peyre

filename = 'nefertiti.off';
[vertex,face] = read_mesh(filename);

options.method = 'freeboundary';
vertex1 = compute_parameterization(vertex,face, options);

A = triangulation2adjacency(face);
plot_graph(A,vertex1);
