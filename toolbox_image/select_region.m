function M1 = select_region(M)

% select_region - extract some sub image
%
%   M1 = select_region(M);
%
%   Copyright (c) 2005 Gabriel Peyré

M1 = M;
[n,p] = size(M);

b = 1;

while b==1
    clf;
    imagesc(M1);
    colormap gray(256);
    axis image; axis off;
    [x,y,b] = ginput(1);
    x = round(x); x = max(min(x,n),1);
    y = round(y); y = max(min(y,p),1);
    M1(y,x) = Inf;
end