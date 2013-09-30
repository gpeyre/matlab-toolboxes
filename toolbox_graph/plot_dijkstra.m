function plot_dijkstra(A, vertex, S, path, start_points,end_points, options )

% plot one state in dijkstra algorithm.
%
%   dijkstra_plot(A, S, vertex, options );
%
%   A is the adjacency matrix of the graph.
%   S is the state after propagation.
%   vertex is the 2D location of the vertices.
%
%   You can add the following options : 
%       'options.far_point_style' : style used for display of non-reached vertices.
%       'options.open_point_style' : style used for display of open vertices.
%       'options.closed_point_style' : style used for display of closed vertices.
%       'options.start_point_style' : style used for display of start vertices.
%       'options.end_point_style' : style used for display of start vertices.
%       'options.path_style' : style used for display of path.
%
%   If 'target' is not empty, it will plot the path back.
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;

if isfield(options, 'far_point_style')
    far_point_style = options.far_point_style;
else
    far_point_style = 'k.';
end

if isfield(options, 'open_point_style')
    open_point_style = options.open_point_style;
else
    open_point_style = 'r.';
end

if isfield(options, 'closed_point_style')
    closed_point_style = options.closed_point_style;
else
    closed_point_style = 'b.';
end

if isfield(options, 'start_point_style')
    start_point_style = options.start_point_style;
else
    start_point_style = 'ro';
end

if isfield(options, 'end_point_style')
    end_point_style = options.end_point_style;
else
    end_point_style = 'go';
end

if isfield(options, 'path_style')
    path_style = options.path_style;
else
    path_style = 'r';
end

if nargin<5
    start_points = [];
end
if nargin<6
	end_points = []
end


if isfield(options, 'graph_style')
    graph_style = options.graph_style;
else
    graph_style = 'k:';
end

if isfield(options, 'point_size')
    point_size = options.point_size;
else
    point_size = 22;
end


O = find(S==0); % open list
C = find(S==-1); % close list

n = size(A,1);

if size(vertex,2)~=2
    vertex = vertex';
end
if size(vertex,2)~=2
    error('Works only for 2D graphs');
end

hold on;
gplot( A,vertex, graph_style );
if ~isempty(far_point_style)
    plot( vertex(:,1), vertex(:,2), far_point_style, 'MarkerSize', point_size );
end
if ~isempty(open_point_style)
    plot( vertex(O,1), vertex(O,2), open_point_style, 'MarkerSize', point_size );
end
if ~isempty(closed_point_style)
    plot( vertex(C,1), vertex(C,2), closed_point_style, 'MarkerSize', point_size );
end
if ~isempty(start_points)
    plot( vertex(start_points,1), vertex(start_points,2), start_point_style );
end
if ~isempty(end_points)
    plot( vertex(end_points,1), vertex(end_points,2), end_point_style );
end

if ~isempty(path)
    if ~iscell(path)
        path = {path};
    end
    for i=1:length(path)
        plot( vertex( path{i}, 1 ), vertex( path{i}, 2 ), path_style )
    end
end

hold off;
axis tight;
axis off;