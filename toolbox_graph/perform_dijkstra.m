function [D,S] = perform_dijkstra(W, start_points, options)

% perform_dijkstra - launch the Dikstra algorithm.
%
%   [D,S] = perform_dijkstra(W, start_points, options);
%
%   If you want a fast dijkstra to compute path from a set of points, use
%       D = perform_dijkstra_fast(sparse(W), point_set);
%   and D(i,j) will be distance between point_set(i) and j.
%
%   'W' is the weight matrix (the highest, the slowest the front will move).
%       W(i,j) is the cost of moving from i to j, 
%       W should be sparse matrix (otherwise 0 entries means no connexion).
%   'start_points' is a 1 x k array, start_points(:,i) is the ith starting
%   point .
%
%   Optional:
%   - You can provide special conditions for stop in options :
%       'options.end_points' : stop when these points are reached
%       'options.nb_iter_max' : stop when a given number of iterations is
%          reached.
%   - You can provide an heuristic in options.H (typically that try to guess the distance
%       that remains from a given node to a given target).
%       This is an array of same size as W.
%
%   Copyright (c) 2005 Gabriel Peyr?

options.null = 0;

if isfield(options,'end_points')
    end_points = options.end_points;
else
    end_points = [];
end

if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 1;
end

if isfield(options,'nb_iter_max')
    nb_iter_max = options.nb_iter_max;
else
    nb_iter_max = Inf;
end

if isfield(options,'H')
    H = options.H;
else
    H = [];
end

if isfield(options,'use_mex')
    use_mex = options.use_mex;
else
    use_mex = 1;
end

nb_iter_max = min(nb_iter_max, 1.2*max(size(W))^2);

W = sparse(W);

% use fast C-coded version if possible
if exist('perform_dijkstra_propagation')~=0 && use_mex
    [D,S] = perform_dijkstra_propagation(W',start_points-1,end_points-1,nb_iter_max, H); % use transposate
else
    [D,S] = perform_dijkstra_propagation_slow(W,start_points,end_points,nb_iter_max, H);
end

% replace C 'Inf' value (1e9) by Matlab Inf value.
I = find( D>1e8 );
D(I) = Inf;
