function D = compute_distance_matrix(X)

% compute_distance_matrix - compute pairwise distance matrix.
%
%   D = compute_distance_matrix(X);
%
%   We have D(i,j)=|X(:,i)-X(:,j)|^2.
%
%   Copyright (c) 2004 Gabriel Peyré

[m,p] = size(X);
X2 = sum(X.^2,1);
D = repmat(X2,p,1)+repmat(X2',1,p)-2*X'*X;