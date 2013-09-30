function T = compute_stick_tf(v,N,n,sigma,c)

% compute_stick_tf - compute the voting kernel of a stick tensor.
%
%   T = compute_stick_tf(v,n,N,sigma,c);
%
%   'v' is the non-null eigenvector of the stick tensor (i.e. the tensor is v*v').
%   'N' is the total size of the image (supposed to lie in [0,1]x[0,1]).
%
%       The equation for the decay of the field is
%           DF(s,k) = exp(-(s^2+c*k^2)/sigma^2)
%       where 's' is the length of the circle that joins
%       the two points, and 'k' is its curvature.
%
%   Optional:
%   'sigma' control the scale of the voting field.
%   'c' control the ratio between distance and curvature.
%       The higher 'c' is, the narrower the field will be.
%   'n' is the size of the kernel (should be odd).
%
%   See also compute_stick_ball.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    error('Not enough arguments');
end
if nargin<3
    n = N/2;
end
if nargin<4
    sigma = 0.2;
end
if nargin<5
    c = 0.1*sigma;
end

% n should be odd
n = floor(n/2)*2+1;
nn = (n-1)/2;
% v should be of unit length
v = v/norm(v,'fro'); v = v(:);
w = [-v(2);v(1)];   % normal vector

x = (-nn:nn)/N;
[Y,X] = meshgrid(x,x);

% rotate the grid to align v to [1;0], ie. multiply by [v,w]'
M = [X(:),Y(:)]';
A = [v,w];
M = A'*M;
X = reshape( M(1,:), [n,n] );
Y = reshape( M(2,:), [n,n] );

% angle
theta = atan2(Y,X);
thetas = abs(theta);
I = find(thetas>pi/2);
thetas(I) = pi-thetas(I);
% length
L = sqrt(X.^2+Y.^2);
% arc-length
s = zeros(n,n);
I = find(L~=0 & thetas~=0);
s(I) = L(I).*thetas(I)./sin(thetas(I));
I = find(L==0 | thetas==0);
s(I) = L(I);
% curvature
k = zeros(n,n);
I = find(L~=0);
k(I) = 2*sin(thetas(I))./L(I);
% attenuation
DF = exp(-(s.^2+c^2*k.^2)/sigma^2);

% the direction vector is a = A*[-cos(2*theta);-sin(theta)]
% and the tensor is DF*a*a'
b = zeros(n,n,2);
b(:,:,1) = -cos(2*theta);
b(:,:,2) = -sin(2*theta);
% rotate point wise
a = b;
a(:,:,1) = A(1,1).*b(:,:,1) + A(1,2).*b(:,:,2);
a(:,:,2) = A(2,1).*b(:,:,1) + A(2,2).*b(:,:,2);
% compute rigidity tensor associated to a
T = zeros(n,n,2,2);
T(:,:,1,1) = DF.*a(:,:,1).^2;
T(:,:,1,2) = DF.*a(:,:,1).*a(:,:,2);
T(:,:,2,1) = T(:,:,1,2);
T(:,:,2,2) = DF.*a(:,:,2).^2;