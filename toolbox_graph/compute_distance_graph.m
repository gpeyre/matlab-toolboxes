function D = compute_distance_graph(W, point_list)

% compute the distance between each point on a graph
%
%   D = compute_distance_graph(W, point_list);
%
%   Uses either 'perform_dijkstra_fast' or 'perform_dijkstra' mex code
%   (depending on which is available).
%
%   D(i,j) is the geodesic graph distance between vertex point_list(i) and
%   vertex j.
%
%   Copyright (c) 2006 Gabriel Peyr?

n = size(W,1); % number of points in the graph
if nargin<2
    point_list = 1:n;
end

%% use the fastest code available
if exist('perform_dijkstra_fast')==3
	D = perform_dijkstra_fast(W, point_list);
    return;
end
    
%% use slow mex code
W = full(W);

D = zeros(length(point_list),n);
hh = waitbar(0,['Computing distances.']);
for i=1:length(point_list)
    waitbar( i/length(point_list) ,hh);
    warning off;
    [d,S] = perform_dijkstra(W, point_list(i));
    warning on;
    D(i,:) = d(:)';
end
close(hh);

% symmetrize
if length(point_list)==n
    D = (D+D')/2;
end