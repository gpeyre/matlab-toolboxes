function h = compute_gabor_filter(n,sigma,theta,f, options)

% compute_gabor_filter - builds a 2D gabor filter.
%
% h = compute_gabor_filter(n,sigma,theta,f);
%
%   theta is the orientation
%   f is the frequency
%   sigma is the width
%
%   Set options.iscomplex=1 if you want a complex Gabor, 
%   otherwise it is real (with cosine modulation).
%
%   If you set theta=[] then the Gabor is radial (no orientation).
%
%   The filter is normalized in norm L1 and with 0 mean.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'iscomplex')
    iscomplex = options.iscomplex;
else
    iscomplex = 0;
end

x = linspace( -n/2,n/2,n );
[Y,X] = meshgrid(x,x);

if not(isempty(theta))
    v = -X*sin(theta) + Y*cos(theta);
else
    v = sqrt(X.^2 + Y.^2);
end

if iscomplex
    Z = exp( 2i*pi*f*v );
else
    Z = cos( 2*pi*f*v );
end

h = Z .* exp(-(X.^2+Y.^2) / (2*sigma^2) );
h = h - mean(h(:));
h = h / sum(abs(h(:)));