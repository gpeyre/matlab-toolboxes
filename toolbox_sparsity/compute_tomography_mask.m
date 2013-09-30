function A = compute_tomography_mask(n,nrays)

% compute_tomography_mask - builds a tomography mask
%
%   A = compute_tomography_mask(n,nrays);
%
%   NB: The original example of Romberg and Candes is with nrays=22.
%
%   Copyright (c) 2008 Gabriel Peyre

Theta = linspace(0,pi,nrays+1); Theta(end) = [];
A = zeros(n);
for theta = Theta
    t = linspace(-1,1,3*n)*n;
    x = round(t.*cos(theta)) + n/2+1; y = round(t.*sin(theta)) + n/2+1;
    I = find(x>0 & x<=n & y>0 & y<=n); x = x(I); y = y(I);
    A(x+(y-1)*n) = 1;
end
