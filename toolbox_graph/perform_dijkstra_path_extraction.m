function path = perform_dijkstra_path_extraction(A,D,end_points)

% perform_dijkstra_path_extraction - extract a shortest path
%
%   path = perform_dijkstra_path_extraction(A,D,end_points);
%
%   A is the adjacency matrix.
%   D is the distance matrix.
%
%   Copyright (c) 2005 Gabriel Peyré

if length(end_points)>1
    path = {};
    for i=1:length(end_points)
        path{i} = perform_dijkstra_path_extraction(A,D,end_points(i));
    end
    return;
end

if D(end_points)==Inf
    warning('end point was not reached');
    I = find(D~=Inf);
    [tmp,j] = min(D(I));
    end_points = I(j);
end

path = [end_points];    % the path
while true
    % select neighbors
    N = find( A(path(end),:)>0 );
    if isempty(N)
        return;
    end
    % find minium distance
    [d,I] = min( D(N) );
    if d>=D(path(end))
        % we are on a minima, stop
        return;
    end
    path(end+1) = N(I);
end