% toolbox_graph - a toolbox for graph manipulation and plotting.
%
%   - A graph vith n vertices is represented by an adjacency matrix A
%   of size nxn where A(i,j)=1 if i is linked to j and A(i,j)=0 otherwise.
%   You can create such a graph using an OFF file (triangulation, see below)
%   or using for instance 
%       [A,vertex] = build_graph_from_image( M, tresh, connect );
%   Note that the combinatorial structure (the matrix A) can
%   be used together with the position of the vertex (the array vertex), 
%   for instance for display use or for geometric computation (flattening, etc).
%   - A triangulation of n vertices and m triangle (faces) is represented
%   via the position of the vertices (which is an array of float of size nx3) 
%   and the connectivity of the faces (which is an array of size mx3).
%   You can create such a triangulation using for example data from an OFF
%   file :
%       [vertex,face] = read_off('venus.off');
%   You can extract the graph structure (i.e. build the adjacency matrix)
%   using :
%       A = triangulation2adjacency(face);
%
%   Here is a list of the relevant functions.
%   Basic manipulation:
%       - compute_vertex_ring - compute the 1 ring of each vertex in a triangulation.
%       - compute_normal - compute the normal at each vertex of a triangulation, using mean of the faces normals.
%       - triangulation2adjacency - compute the adjacency matrix of a given triangulation.
%       - reverse_orientation - reverse the orientation of the face of a mesh.
%   Graph and triangulation creation:
%       - gen_cyclic_graph - generate a Caley graph associated with generating set S in Z/nZ.
%       - gen_square_graph - generate a Caley graph on a square grid.
%       - build_graph_from_image - build a graph from an image.
%       - compute_base_mesh - generate a simple triangulation.
%       - read_off - read data from OFF file.
%       - perform_mesh_subdivision - perform a 1:4 subdivision.
%       - subdivide_sphere - subdivide each triangle into 4 smaller triangles, and then project on a sphere.
%   Parameterization and flattening methods:
%       - dijkstra - find shortest paths in graphs.
%       - floyd - compute shortest distances on a graph using Floyd algorithm.
%       - build_euclidean_weight_matrix - build a weighted adjacenty matrix that reflect position of the vertices.
%       - isomap - computes Isomap embedding using the algorithm of Tenenbaum, de Silva, and Langford (2000). 
%       - compute_mesh_laplacian - return a laplacian of a given triangulation (can be combinatorial or geometric).
%       - compute_parametrization - compute a planar parameterization of a given disk-like triangulated manifold.
%		- compute_mesh_gradient, compute_mesh_weights - compute gradient and low pass operator on meshes.
%   Display
%       - off_previewer - load a file in OFF/PLY/WRL/SMF format and display the model.
%       - plot_mesh - plot a 3D mesh.
%
%   Copyright (c) 2004 Gabriel Peyré

disp('type ''help toolbox_graph''.');