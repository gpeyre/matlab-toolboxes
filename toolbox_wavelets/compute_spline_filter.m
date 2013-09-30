function [g,h] = compute_spline_filter(s,N)

% compute_spline_filter - compute the spline filter of a wavelete transform
%
% function [g,h] = compute_spline_filter(s,N);
%
%    'N' is the length of the filter (must be even).
%    's' is the order of the spline.
%    'g' is the low pass filter.
%    'h' is the high pass filter.
%


x = (0:N-1)*2*pi/(N-1);
gg = ; % compute here the fourier transform
g = ifft(gg);
g = [g(N/2+1:end),g(1:N/2)];
h = (-1).^(0:N-1) .* g;