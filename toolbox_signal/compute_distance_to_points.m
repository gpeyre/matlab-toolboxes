function D = compute_distance_to_points(X,Y, options)

% compute_distance_to_points - compute euclidean distance to a set of points.
%
%   D = compute_distance_to_points(X,Y)
%
%   'X' is a [d,n] matrix, X(:,i) is the ith point living in R^d.
%   'Y' is a [d,k] matrix (in practice, k<<n).
%   D(i,j) = |X(:,j)-Y(:,i)|^2.
%
%   Copyright (c) 2004 Gabriel Peyré

m = size(Y,2);
n = size(X,2);
D = zeros(m,n);
d = size(X,1);

% dimension reduction
options.null = 0;
dr = getoptions(options, 'dr', d);
if dr<d
    nbexemplars = min(n,5000);
    sel = randperm(n);
    sel = sel(1:nbexemplars);
    % compute PCA
    [P,X1,v,Psi] = pca(X(:,sel),dr);
    % perfrom PCA projection
    X = X - repmat( Psi, [1 n] ); X = P'*X;
    Y = Y - repmat( Psi, [1 m] ); Y = P'*Y;
end

for k=1:m
    % distance to seed
    D(k,:) = sum( (X - repmat(Y(:,k),1,n)).^2 );
end