function [D,nn_list] = compute_nn_distance(X,nbr_nn, options)

% compute_nn_distance - compute the distance to the nearest neighbors
%
%   [D,nn_list] = compute_nn_distance(X,nbr_nn, options);
%
% X is a (d,n) set of n points in R^d
% nbr_nn is the number of nearest neighbors to retrieve.
%
%   You can set Y=options.Y, otherwise the algorithm assumes
%       the target points are Y=X.
%
% D is a (n,nbr_nn) matrix of distance, D(i,j) is distance between point
%   Y(:,i) and point X(nn_list(i,j),:)
% nn_list(i,:) is the set of nearst neighbors
%
%   options.use_nntools = 0 force the use of the slow matlab code
%       for nearest neighbors computations.
%
%   options.exlude_self = 1 to avoid taking a point it self neighbor
%       (works only when Y=X).
%
%   If you set options.pca_numvecs < d, then a preprocessing step of
%       dimensionnality reduction is performed
%
%   Copyright (c) 2006 Gabriel Peyr?

options.null = 0;

if isfield(options, 'use_nntools')
    use_nntools = options.use_nntools;
else
    use_nntools = 1;
end
if isfield(options, 'exlude_self')
    exlude_self = options.exlude_self;
else
    exlude_self = 0;
end

if isfield(options, 'Y')
    Y = options.Y;
else
    Y = [];
end
if exlude_self && ~isempty(Y)
    warning('You can not use exlude_self with options.Y enabled.');
end

if nargin<2
    nbr_nn = size(X,2);
end


nbr_nn = min(nbr_nn,size(X,2));

d = size(X,1);
n = size(X,2);
if d>2*n
    warning('Matrix seems to be of wrong dimension, should be dimension x nbr_points');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensionnality reduction for fast search
if isfield(options, 'pca_numvecs')
    m = options.pca_numvecs;
else
    m = d;
end
if m>0 && m<d
    %%% perform dimensionnality reduction 
    [P,tmp,v,Psi] = pca(X,m);
    X = X - repmat(Psi,[1,size(X,2)]);
    X = P' * X;
    if ~isempty(Y)
        Y = Y - repmat(Psi,[1,size(Y,2)]);
        Y = P' * Y;
    end
    d = m;
end

if not(exist('nn_prepare')==3) && use_nntools
    warning('TSTool not installed, use slow Matlab code instead.');
    use_nntools = 0;
end


if isempty(Y)
    %%% compute distance from X to X %%%
    if exlude_self==0
        nbr_nn = nbr_nn - 1;
    end
    if use_nntools
        % use fast mex code
        atria = nn_prepare(X');
        [nn_list,D] = nn_search(X', atria, 1:n, nbr_nn, 0);
    else
        % use slow matlab code
        D1 = sqrt( compute_distance_matrix(X) );
        D1 = D1 + diag( Inf + zeros(size(D1,1),1) );
        % find closest points
        [D,nn_list] = sort(D1);
        D = D(1:nbr_nn,:)';
        nn_list = nn_list(1:nbr_nn,:)';
    end
    if exlude_self==0
        % add self reference
        nn_list = [(1:n)', nn_list];
        D = [zeros(n,1), D];
    end
else
    %%% compute distances from Y to X %%%
    if exist('nn_prepare')>0 && use_nntools
        % use fast mex code
        atria = nn_prepare(X');
        [nn_list,D] = nn_search(X', atria, Y', nbr_nn, 0);
    else
        % use slow matlab code
        D1 = sqrt( compute_distance_matrix(X,Y) );
        % find closest points
        [D,nn_list] = sort(D1,2);
        D = D(:,1:nbr_nn)';
        nn_list = nn_list(:,1:nbr_nn)';
    end
end