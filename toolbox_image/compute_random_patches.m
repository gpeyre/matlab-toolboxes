function H = compute_random_patches(M,w,m, wmax)

% compute_random_patches - extract patches
%
%   H = compute_random_patches(M,w,m, wmax);
%
%   w is the width of the patches, m the number of patches.
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin<4
    wmax = w;
end

sigma = w/wmax;
if mod(sigma,1)~=0
    error('wmax should be a multiple of w');
end

n = size(M,1);

if sigma>1
    h = compute_gaussian_filter([17 17],sigma/(2*n),[n n]);
    M = perform_convolution(M,h);
end

X = floor(rand(m,1)*(n-w-1))+1;
Y = floor(rand(m,1)*(n-w-1))+1;

X = reshape(X, [1 1 1 m]);
Y = reshape(Y, [1 1 1 m]);


% sampling locations
[dY,dX] = meshgrid(0:w-1,0:w-1);
IX = repmat(dX,[1 1 1 m]) + repmat( X, [w w 1 1] );
IY = repmat(dY,[1 1 1 m]) + repmat( Y, [w w 1 1] );
IX = mod(IX-1,n)+1;
IY = mod(IY-1,n)+1;
I = sub2ind([n n], IX, IY);


for s=1:size(M,3)
    U = M(:,:,s);
    H(:,:,s,:) = U(I);
end

if sigma>1
    % sub-samples
    H = H(1:sigma:end,1:sigma:end,:,:);
end