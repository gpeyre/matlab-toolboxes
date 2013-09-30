function W = build_euclidean_weight_matrix(A,vertex,pad_with)

% build_euclidean_weight_matrix - build a weighted adjacenty matrix
%   that reflect position of the vertices.
%   The weight is Inf if vertices are not connected.
%   You can put something different from Inf (e.g. 0) 'using pad_with'.
%
%   W = build_euclidean_weight_matrix(A,vertex [,pad_with] );
%
%   Copyright (c) 2004 Gabriel Peyr?

if nargin<2
    error('Not enough arguments.');
end
if nargin<3
    pad_with = Inf;
end

if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

n = size(A,1);
if size(vertex,1)<n
    error('Not enough vertices');
end

W = sparse(n,n);
W(:,:) = pad_with;

for i = 1:n
    for j=1:n
        if A(i,j)~=0
            W(i,j) = norm( vertex(i,:)-vertex(j,:), 'fro' );
        end
    end    
end