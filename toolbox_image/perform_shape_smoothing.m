function M = perform_shape_smoothing(M,sigma)

% perform_shape_smoothing - perform morphological smoothing
%
% M = perform_shape_smoothing(M,sigma);
%
%   M is a binary shape (1 is assumed to the shape and 0 the background).
%   sigma is the width of the smoothing, in pixels.
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin<2
    sigma = 2;
end

n = size(M,1);

M = double(rescale(M)>0.5);
% smooth
h = compute_gaussian_filter( round([n n]/8)*2+1,sigma/(2*n),[n n]);
M = perform_convolution(M,h);
% threshold
M = double(rescale(M)>0.5);