function xy = perform_geodesic_embedding(X,options)

% perform_geodesic_embedding - compute embedding using distance to seed points.
%
%   xy = perform_geodesic_embedding(X,options);
%
%   You can provide:
%       options.landmarks: the list of the numbers of the 
%           extremal points used to compute the embedding.
%           The dimension of the embedding is m=(nbr.landmarks)/2
%       options.dim: the embedding dimension. The dim*2 needed lamdmarks
%           are then computed automatically using farthest points.
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;

if isfield(options, 'landmarks')
    landmarks = options.landmarks;
    m = length(landmarks)/2;
else
   if isfield(options, 'dim') 
       m = options.dim;
   else
       m = 2;
   end
end

if isempty(landmarks)
    % compute automatically extremal points
    error('Not yet implemented');
end

ld1 = landmarks(1:2:end);
ld2 = landmarks(2:2:end);
n = size(X,2);

options.nn_nbr = 8;
G = compute_nn_graph(X,options);

Ws = G;
Ws(find(Ws==Inf)) = 0;
Ws = sparse(Ws);
D = dijkstra_fast(Ws, landmarks);

xy = zeros(m,n);

for i=1:m
    i1 = 2*i-1;
    i2 = 2*i;
    xy(i,:) = D(i1,:) ./ ( D(i1,:) + D(i2,:) );    
end