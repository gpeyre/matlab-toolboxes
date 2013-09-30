function D = dist_eucl(x,y)

% dist_eucl - compute the euclidean distances between points.
%
%   D = dist_eucl(x,y);
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(x,2);
p = size(y,2);
d = size(x,1);

[J,I] = meshgrid(1:p,1:n);

D = zeros(1,n*p);
for s=1:d
    D = D + (x(s,I)-y(s,J)).^2;
end
D = sqrt(D);
D = reshape(D,n,p);