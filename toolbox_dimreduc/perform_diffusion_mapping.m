function [Y,s] = perform_diffusion_mapping(X, dim, options)

% perform_diffusion_mapping - perform the diffusion dimension reduction
%
%   [Y,s] = perform_diffusion_mapping(X, dim, options);
%
%   X = data as D x N matrix (D = dimensionality, N = #points)
%   dmax = max embedding dimensionality
%   Y = embedding as dim x N matrix
%   s = embedding spectra.
%
%   You can compute a diffusion embedding for varying time step t using
%       Yt = Y .* repmat( s.^t, [N 1] )
%
%   Compute the embedding as the first eigenvectors of the
%   symmetric transition matrix H defined by the following procedure:
%
%   1) Kernel K is defined by
%       K_ij = exp( -|X(:,i)-Y(:,j)|^2 / sigma^2 )
%       where sigma^2 = 1/N \sum_i min_{j != i} |X(:,i)-Y(:,j)|^2
%   2) If options.normalize==1, then the kernel is normalized 
%       p = K * ones(N,1)
%       K = k ./ (p*p')
%   3) The matrix is made approximatrly stochastic and symmetric
%       v = sqrt( K * one(N,1) )
%       H = K ./ (v*v')
%
%   The algorithm is described in 
%       R.R.Coifman, S. Lafon, A.B. Lee, M. Maggioni, B. Nadler, F. Warner and S.W. Zucker, 
%       “Geometric diffusions as a tool for harmonic analysis and structure definition of data”, 
%       Proceedings of the National Academy of Sciences, Vol 102(21), May 2005
%
%   Copyright (c) 2006 Gabriel Peyré


options.null = 0;
if nargin<2
    dim = 2;
end

N = size(X,2);

% compute NN relationship
nbr_nn = 20;
[D,nn_list] = compute_nn_distance(X,nbr_nn);

% compute gaussian kernel
sigma = mean( mean(D) )^2;
sigma = 4 * 0.1 * mmax(X)/5;

K = zeros(N);
I = sub2ind([N N], repmat((1:N)', [1 nbr_nn]), nn_list );
K(I) = exp( -D.^2 / sigma );
K = (K+K')/2;
K = sparse(K);

% perform normalization
if isfield(options, 'normalize') && options.normalize==1
    p = K * ones(N,1);
    K = K ./ (p*p');
end
v = sqrt( K * ones(N,1) );
K = K ./ (v*v');

% perform spectral factorization
[U,S,V] = svds(K, dim+1);
U = U ./ repmat(U(:,1), [1 dim+1]);
Y = U(:,2:dim+1);
s = diag(S); s = s(2:end); s = s(:)';