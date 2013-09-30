function M = callback_blurring( M, dir, options )

% callback_blurring - callback for the deconvolution with wavelets
%
%   M = callback_blurring( M, dir, options );
%
%   Copyright (c) 2007 Gabriel Peyre

if isfield(options, 'eta')
    eta = options.eta;
else
    eta = 5;
end
if isfield(options, 'eta')
    Jmin = options.Jmin;
else
    Jmin = 4;
end

n = sqrt(size(M,1));
M = reshape(M,n,n);

if dir==1
    % wavelet synthesis and then blurring
    M = perform_wavelet_transform(M, Jmin, -1, options );
    M = perform_blurring( M, eta );
else
    % transpose
    M = perform_blurring( M, eta );
    M = perform_wavelet_transform(M, Jmin, 1, options );
end

M = M(:);