function xy = Isomap(D,ndims,options); 

% isomap - computes Isomap embedding using the algorithm of 
%             Tenenbaum, de Silva, and Langford (2000). 
%
%   xy = isomap(X,ndims,options); 
%       or 
%   xy = isomap(D,ndims,options); 
%
%   'X' is a D x N matrix (D = dimensionality, N = #points)
%   'D' is either a local distance matrix (with Inf when no connection),
%       or a global one (i.e. it contains already computed geodesic
%       distances between pair of points).
%   OPTIONAL:
%   'ndims' is the number of output dimensions (default =2).
%   'options.distance_mode' is either 'local' (set if D is a local distance matrix)
%       or 'global' (set if D is a geodesic distance).
%   'options.nn_epsilon' : set it if you want to compute the nearest
%       neighbor by thresholding.
%   'options.nn_nbr' : set it if you want to compute the nearest
%       neighbor by using a fixed number of neighbors.
%       
%   Modified by Gabriel Peyr?? from the original code :
%
%    BEGIN COPYRIGHT NOTICE
%
%    Isomap code -- (c) 1998-2000 Josh Tenenbaum
%
%    This code is provided as is, with no guarantees except that 
%    bugs are almost surely present.  Published reports of research 
%    using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%      J. B. Tenenbaum, V. de Silva, J. C. Langford (2000).  A global
%      geometric framework for nonlinear dimensionality reduction.  
%      Science 290 (5500): 2319-2323, 22 December 2000.  
%
%    Comments and bug reports are welcome.  Email to jbt@psych.stanford.edu. 
%    I would also appreciate hearing about how you used this code, 
%    improvements that you have made to it, or translations into other
%    languages.    
%
%    You are free to modify, extend or distribute this code, as long 
%    as this copyright notice is included whole and unchanged.  
%
%    END COPYRIGHT NOTICE

options.null = 0;

if nargin<2
    ndims = 2;
end

if isfield(options, 'distance_mode')
    distance_mode = 'local';
else
    if ~isempty(find(isinf(D)));
        distance_mode = 'local';
    else
        distance_mode = 'global';
    end
end


if isfield(options, 'verb')
    verb = options.verb;
else
    verb = 1;
end


use_landmarks = 0;
points_list = 1:size(D,2);
if isfield(options, 'landmarks')
    use_landmarks = 1;
    landmarks = options.landmarks;
    nbr_landmarks = length(landmarks);
    points_list = landmarks;
end

if size(D,1)~=size(D,2)
    X = D;
    % D is the location of point in space, compute the distance matrix
    % using some NN-graph
    D = compute_nn_graph(X,options);
    distance_mode = 'local';
end

N = size(D,1);

if strcmp(distance_mode, 'local')
    %%%%% Compute shortest paths %%%%%
    if verb
        disp('- computing shortest paths');
    end
    % D = floyd(D);
    Ws = D;
    Ws(Ws==Inf) = 0;
    Ws = sparse(Ws);
    % compute the distance between point_list and the remaining points
    D = compute_distance_graph(Ws, points_list);
    if use_landmarks
        Dfull = D;
        D = D(:,landmarks);
        Nfull = N;
        N = nbr_landmarks;
    end
end

if sum(D==Inf)>0
    warning('Distance matrix contains Inf value (non-connected graph)');
    D(D==Inf) = max(D(not(D==Inf)))*2;
end

%%%%% Construct low-dimensional embeddings (Classical MDS) %%%%%
if verb
    disp('- constructing low-dimensional embeddings.'); 
end
opt.disp = 0; opt.isreal = 1; opt.issym = 1; 
M = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));
[xy, val] = eigs(M, ndims, 'LR', opt); 

for i=1:ndims
    xy(:,i) = xy(:,i)*sqrt(val(i,i));
end

if use_landmarks
    % interpolation on the full set of points
    % x = 1/2 * (L^T) * ( delta_n-delta_x )
    xy1 = zeros(Nfull,ndims);
    % transpose of embedding
    LT = xy'; 
    for i=1:ndims
        LT(i,:) = LT(i,:) / val(i,i);
    end
    deltan = mean(D,2);
    for x=1:Nfull
        deltax = Dfull(:,x);
        xy1(x,:) = 1/2 * ( LT * ( deltan-deltax ) )';
    end
    xy = xy1;
end