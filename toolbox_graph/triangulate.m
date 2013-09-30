% triangulate - interface to 'triangle' library. 
%
%   [vertex, face] = compute_constrained_delaunay(vertex, edges, holes);
%
%   Triangulates a set of points with boundary, using a constrainted Delaunay 
%   triangulation. It uses the Triangle library of Jonathan Richard Shewchuk:
%       http://www-2.cs.cmu.edu/~quake/triangle.html
%
%   Copyright (c) 2004 Gabriel Peyré

function [nodes, triangles] = triangulate(nodes, boundary, holes)

if size(nodes,1)<size(nodes,2)
    nodes = nodes';
end
if size(boundary,1)<size(boundary,2)
    boundary = boundary';
end
if size(holes,1)<size(holes,2)
    holes = holes';
end

nodes = [(1:size(nodes,1))', nodes];
boundary = [(1:size(boundary,1))', boundary];
holes = [(1:size(holes,1))', holes];

n = size(nodes, 1);
m = size(boundary, 1);
h = size(holes, 1);

% Write out domain in PSLG format.
fid = fopen('mesh.poly', 'w');
if fid<0
    error('Unable to create file mesh.poly.');
end
fprintf(fid, '%d 2 0 0\n', n);
for i = 1:n
    fprintf(fid, '%d %f %f\n', nodes(i,1), nodes(i,2), nodes(i,3));
end
fprintf(fid, '%d 0\n', m);
for j = 1:m
    fprintf(fid, '%d %d %d\n', boundary(j,1), ...
        boundary(j,2), ...
        boundary(j,3));
end
fprintf(fid, '%d\n', h);
for k = 1:h
    fprintf(fid, '%f %f %f\n', holes(k,1), holes(k,2), holes(k,3));
end
fclose(fid);

% Triangulate domain.
!triangle -pq mesh.poly > /dev/tmp

% Read in new nodes.
fid = fopen('mesh.1.node', 'r');
if fid<0
    error('No file mesh.1.node produced.');
end
n = fscanf(fid, '%d', 1);
dummy = fscanf(fid, '%d', 3);
nodes = [];
for i = 1:n
    point = [ fscanf(fid, '%d', 1)' fscanf(fid, '%f', 2)' ];
    nodes = [ nodes; point ];
    fscanf(fid, '%d', 1);
end
fclose(fid);

% Read in triangle mesh.
fid = fopen('mesh.1.ele', 'r');
m = fscanf(fid, '%d', 1);
dummy = fscanf(fid, '%d', 2);
triangles = [];
for j = 1:m
    triangles = [ triangles; fscanf(fid, '%d', 4)' ];
end
fclose(fid);

triangles = triangles(:,2:end)';
nodes = nodes(:,2:end)';

% Delete temporary files.
!rm -f mesh.poly mesh.1.poly mesh.1.node mesh.1.ele
