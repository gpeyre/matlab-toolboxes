function [LocP,LocPsi,Segm] = perform_local_pca(X, options)

% perform_local_pca - perform local principal component
%
% [LocP,LocPsi,Segm] = perform_local_pca(X, options);
%
%   X is the data set of size (m,p), X(:,i) is the ith point 
%       leaving in dimension p.
%
%   LocP(:,:,i) is the matrix whose d columns are the basis
%       of the tangent plane in the ith cluster.
%   LocPsi(:,i) is m dimensional vector which is the center 
%       of the ith cluster.
%   Segm(i) tels to which cluster belongs the ith point of X. 
%
%   You can define
%       options.nbr_clusters : number of center points
%       options.nb_iter : number of iteration for the clustering.
%       options.nbr_neighbors : number of NN for tangent planes estimation
%       options.dim=d : dimensionality of tangent planes
%
%   Only a k-means clustering is used to perform segmentation.
%   
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;
if isfield(options, 'nbr_clusters')
    p1 = options.nbr_clusters;
else
    p1 = 200;
end
if ~isfield(options, 'nb_iter')
    options.nb_iter = 4;
end
if isfield(options, 'nbr_neighbors')
    k = options.nbr_neighbors;
else
    k = 100;
end
if isfield(options, 'dim')
    d = options.dim;
else
    d = 6;
end

m = size(X,1);
p = size(X,2);
q = 1000; % number of trial point to estimate centers
fprintf('Performing k-mean ... ');
Ind = floor(rand(q,1)*p)+1;
[B,Ys,E] = perform_kmeans(X(:,Ind),p1,options);
fprintf(['done, energy=' num2str(E,3) '.\n']);

% compute manifold segmentation
atria = nn_prepare(Ys');
[Segm,distance] = nn_search(Ys', atria, X', 1, 0);

% estimate some tangent plane
fprintf('--> Compute NN ... ');
atria = nn_prepare(X');
[Neigh,distance] = nn_search(X', atria, Ys', k, 0);
fprintf('done.\n');

% compute local pca
fprintf('--> Compute local PCA ... ');
LocP = zeros( m,d,p1 ); % local projection matrix
LocPsi = zeros( m,p1 ); % local mean
for i=1:p1
    ind = Neigh(i,:); % the k NN
    [LocP(:,:,i),tmp,v,LocPsi(:,i)] = pca(X(:,ind),d);
end
fprintf('done.\n');