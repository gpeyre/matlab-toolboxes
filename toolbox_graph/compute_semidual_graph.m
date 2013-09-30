function [Ad,vertexa] = compute_semidual_graph(face, vertex)

% compute_semidual_graph - compute a half dual graph
%
%   [A,vertex1] = compute_semidual_graph(face, vertex);
%
%   face is the face structure of the triangulation.
%   vertex is optional and gives location of the vertices
%   vertex1 is the position of the new vertex.
%   A is the adjacency matrix of the graph
%
%   This function creates a graph that links the center of the faces
%   of the triangulation to the vertex of the triangulation.
%   The new number of vertices is nverts+nfaces and the new
%   number of faces is nfaces*3.
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin<2
    vertex = [];
end
% the code assume that face is of size (nface,3)
[vertex,face] = check_face_vertex(vertex,face);

nfaces = size(face, 2);
nverts = max(face(:));

% generate a  graph that links center of faces to vertices
if nargin>1
    vertexc = [   ...
        sum(reshape(vertex(1,face),[3 nfaces]), 1)/3; ...
        sum(reshape(vertex(2,face),[3 nfaces]), 1)/3; ...
        sum(reshape(vertex(3,face),[3 nfaces]), 1)/3 ];
else
    vertexc = [];
end

% both primal and dual vertices
vertexa = [vertex, vertexc];
nvertsa = nverts+nfaces;
Ad = zeros(nvertsa);
for i=1:nfaces
    iface = i+nverts; % index of the center of the faces
    Ad(iface,face(:,i)) = 1;
end
Ad = max(Ad,Ad'); % symetrize