function D = dist_sphere(x,y)

% dist_eucl - compute the spherical distances between points.
%
%   D = dist_sphere(x,y);
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(x,2);
p = size(y,2);
d = size(x,1);

% renormalize
d = sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2);
x(1,:) = x(1,:)./d;
x(2,:) = x(2,:)./d;
x(3,:) = x(3,:)./d;

d = sqrt(y(1,:).^2+y(2,:).^2+y(3,:).^2);
y(1,:) = y(1,:)./d;
y(2,:) = y(2,:)./d;
y(3,:) = y(3,:)./d;

[J,I] = meshgrid(1:p,1:n);

D = x(1,I).*y(1,J) + x(2,I).*y(2,J) + x(3,I).*y(3,J);

D = real( acos(D) ) ;
D = reshape(D,n,p);