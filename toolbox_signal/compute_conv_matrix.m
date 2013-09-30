function y = compute_conv_matrix(x)

% compute_conv_matrix - compute a circular convolution matrix.
%
% M = compute_conv_matrix(x);
%  
%     Copyright (c) 2005 Gabriel Peyré

x = x(:);
n = length(x);
[Y,X] = meshgrid(0:n-1,0:n-1);
Z = mod( X+Y, n )+1;
y = x(Z);