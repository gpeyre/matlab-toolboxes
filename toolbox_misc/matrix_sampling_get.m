function y = matrix_sampling_get(M, I)

% matrix_sampling_get - get the value of matrix M at given index.
%
%   y = matrix_sampling_get(M, I);
%
%   M is a 2D array, and I=[I1;I2] is of size 2xn. Return
%   y(i) = M( I1(i), I2(i) )
%
%   Copyright (c) 2004 Gabriel Peyré

J = sub2ind(size(M), I(1,:),I(2,:) );
y = M(J);
