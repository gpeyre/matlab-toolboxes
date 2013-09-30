function M = perform_median_filtering(M,k)

% perform_median_filtering - perform moving average median
%
%   M = perform_median_filtering(M,k);
%
%   k is the half width of the window (detult k=1).
%
%   This filtering method is very efficient to remove impulsive or salt and
%   peper noise.
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<2
    k = 1;
end
w = 2*k+1;
n = size(M,1);

options.sampling = 'uniform';
H = compute_patch_library(M,k,options);
H = reshape(H, [w*w n*n]);
H = median(H); 
M = reshape(H, n,n);