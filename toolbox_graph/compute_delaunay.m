function face = compute_delaunay(vertex)

% compute_delaunay - compute faces of a Delaunay triangulation
%
%   face = compute_delaunay(vertex);
%
%   Copyright (c) 2008 Gabriel Peyre

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end

face = delaunay(vertex(1,:), vertex(2,:))';