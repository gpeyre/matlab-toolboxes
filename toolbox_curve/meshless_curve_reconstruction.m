function [I,t] = meshless_curve_reconstruction(c, phi,R)

% meshless_curve_reconstruction - reconstruct a curve from point cloud.
%
%   [I,t] = meshless_curve_reconstruction(c, phi, R);
%
%   I is the ordering of the points in c (a 2 x npts matrix).
%   t is the parameterization of the set of points, t(i) \in [0,1]
%   is the parameter for point c(:,i).
%   phi is a callback function eg. phi(x)=1/x if x<1 and 0 otherwise.
%   R is the locallity radius (the smaller, the sparser).
%
%   Copyright (c) 2004 Gabriel Peyré

npts = size(c,2);

x = c(1,:);
y = c(2,:);
if nargin<3
    R = abs(max(x)-min(y)) + abs(max(x)-min(y));
    R = R*10/npts;
end
if nargin<2
    phi = inline('1./t .* (t<1) .* (t>0)');
end

% select a start point and an end point
[tmp,istart] = min(x);
[tmp,iend] = max(x);

A = zeros(npts,npts);
b = zeros(npts,1);
for i=1:npts
    if i==istart
        A(i,i) = 1;
        b(i) = 0;
    elseif i==iend
        A(i,i) = 1;
        b(i) = 1;
    else
        % distance to current point
        d = sqrt( ( x-x(i) ).^2 + ( y-y(i) ).^2 );
        d(i) = Inf;
        d = phi(d/R);        
        A(i,:) = -d;
        A(i,i) = sum( d );
        if 0
            I1 = find(d>0);
            I2 = find(d==0);
            clf;
            hold on;
            plot(x(I1), y(I1), 'r.')
            plot(x(I2), y(I2), 'r.');
            hold off;
        end
    end
end

% solve system
t = A\b;
[tmp,I] = sort(t);
t = rescale(t);