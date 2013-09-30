function M1 = perform_image_rotation(M,theta)

% perform_image_rotation - rotate the image
%
%   M1 = perform_image_rotation(M,theta);
%
%   Rotate anti-clockwise.
%
%   Copyright (c) 2008 Gabriel Peyre

n = size(M,1);
x = linspace(-1,1,n);
[X1,Y1] = meshgrid(x,x);
X = cos(theta)*X1 + -sin(theta)*Y1;
Y = sin(theta)*X1 + cos(theta)*Y1;
M1 = interp2(X1,Y1,M,X,Y);
M1(isnan(M1)) = 0;