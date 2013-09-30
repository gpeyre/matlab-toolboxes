function M = perform_fast_rbf_interpolation(v,points, n)

% perform_fast_rbf_interpolation - interpolate using quick and dirty RBF.
%
% M = perform_fast_rbf_interpolation(v, points, n);
%
%   points(:,i) is the ith point with coordinates in 1...n.
%   v(i) is the value of the function.
%
%   This code uses 1/x^3 RBF and force interpolation.
%
%	Works only on images (2D data)
%
%   Copyright (c) Gabriel Peyré 2007

v = v(:);
m = length(v);
D = zeros(n,n,m);
[Y,X] = meshgrid(1:n,1:n);
% compute distance
for i=1:m
    D(:,:,i) = sqrt( (points(1,i)-X).^2 + (points(2,i)-Y).^2 );
end
eta = 0.1;
D = 1./(eta^3+D.^3);
D = D ./ repmat( sum(D,3), [1 1 m] );
% interpolate
M = repmat( reshape(v,[1 1 m]), [n n 1] ) .* D;
M = sum(M,3);