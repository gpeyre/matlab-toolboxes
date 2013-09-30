function [H,P,Psi] = perform_lowdim_embedding(M,options)

% perform_lowdim_embedding - perform a patch wise dimension extension
%
%   [H,options.P, options.Psi] = perform_lowdim_embedding(M,options);
%
%   M = perform_lowdim_embedding(H,options);
%
%   This function lift each pixel of an image (can be a color image)
%   to a vector that incorporate neighboorhood relationship.
%
%   Each pixel is replaced by the vector containing the values of the
%   neighbooring pixels and then dimension reduction is applyed to 
%   avoid manipulating very high dimensional vectors.
%
%   options.ndims gives the dimensionality for PCA.
%
%   Copyright (c) 2006 Gabriel Peyr?


dir = 1;
if size(M,3)>5
    dir = -1;
end

options.null = 0;

k = getoptions(options, 'k', 2);
mask = getoptions(options, 'mask', 'cst');

switch mask
    case 'cst'
        phi = ones(2*k+1);
    case 'linear'
        x = 1-abs(linspace(-1,1,2*k+3)); x = x(2:end-1);
        phi = x'*x;
end

[m,n,s] = size(M);

if dir==1
    ndims = getoptions(options, 'ndims', 25);
    % perform patch wise embedding
    s = size(M,3);
    % extract patches
    options.sampling = 'uniform';
    H = compute_patch_library(M,k,options);
    H = H .* repmat( phi, [1 1 s n*m] );
    % turn into collection of vectors
    H = reshape(H, [s*(2*k+1)^2 n*m]);
end

ndims = min(ndims,size(H,1));

if dir==-1 
    if not(isfield(options, 'P')) || not(isfield(options, 'Psi'))
        error('You must provide options.P and options.Psi.');
    end
    P = options.P; Psi = options.Psi;
    ndims = size(P,2);
    % do un-projection
    M = reshape(M, [n^2 ndims]);
    M = shiftdim(M,1);
    M = (P*M) + repmat( Psi, [1 n^2] );
    % sample at center
    s = size(M,1)/(2*k+1)^2; % number of colors
    M = reshape( M, [s (2*k+1) (2*k+1) n n] );
    M = M(:,(end+1)/2,(end+1)/2,:,:);
    H = squeeze(M);
elseif dir==1
    % compute PCA projection
    if not(isfield(options, 'P')) || not(isfield(options, 'Psi')) || isempty(options.P) || isempty(options.Psi)
        nbexemplars = min(n*m,5000);
        sel = randperm(n*m);
        sel = sel(1:nbexemplars);
        [P,X1,v,Psi] = pca(H(:,sel),ndims);
    else
        P = options.P; Psi = options.Psi;
    end
    ndims = size(P,2);
    % perform actual PCA projection
    H = H - repmat( Psi, [1 n*m] );
    H = P'*H;
    % reshape matrix
    H = reshape(H, [ndims m n]);
    H = shiftdim(H,1);
end
