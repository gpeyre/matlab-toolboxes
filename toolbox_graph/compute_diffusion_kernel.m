function K = compute_diffusion_kernel(A,normalization)

% compute_normalized_kernel - compute a diffusion kernel from a weight/adjacency matrix
%
%   K = compute_diffusion_kernel(A);
%
%   A is a (n,n) adjacency matrix or a weight matrix
%       A(i,j) is the weigth between points, e.g.
%       A(i,j) = exp( -|xi-xj|^2/sigma^2 )
%   if normalization==1, perform full normalization (ie. K will be symmetric)
%   otherwise, just do row equalisation (ie. the row will sum to one, the 
%   matrix is stochastic)
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    normalization = 1;
end

N = size(A,1);

if normalization==1
    D = sparse([1:N], [1:N], 1./sum(A), N, N, N);

    A = D * A * D;

    D = sum(A);
    for j=1:N
        for k=1:j
            A(j,k) = (1/sqrt(D(j))) * (1/sqrt(D(k))) * A(j,k);
            A(k,j) = A(j,k);
        end
    end
else
    D = sparse(diag(1./sum(A)));
    A = D * A;
end

K = sparse(A);