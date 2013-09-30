function H = compute_directional_kernel(sigma1,sigma2,theta,nmax)

% compute_directional_kernel - compute oriented gaussian.
%
% h = compute_directional_kernel(s1,s2,t,nmax);
%
%   s1 is variance along first axis in pixels.
%   s2 is variance along secong axis.
%   t is the angle.
%   nmax is the maxium size of the kernel (optional).
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin<4
    nmax = 128;
end

% size of the kernel
s = max([sigma1(:);sigma2(:)]);
m = min(s*4,nmax);
% m = nmax;
m = round(m/2)*2+1;

p = length(sigma1(:)); % number of filters
H = zeros(m,m,p);

for i=1:p

    % renormalize so that it correspond to pixel size
    s1 = sigma1(i)/2;
    s2 = sigma2(i)/2;
    t = theta(i);

    x = linspace(-m/2,m/2,m);
    [Y,X] = meshgrid(x,x);
    % rotate
    X1 = cos(t)*X + sin(t)*Y;
    Y1 = - sin(t)*X + cos(t)*Y;
    X1 = X1 / s1;
    Y1 = Y1 / s2;
    h = exp( -(X1.^2 + Y1.^2) );
    h = h/sum(h(:));
    H(:,:,i) = h;
end