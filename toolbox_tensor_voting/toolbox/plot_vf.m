function plot_vf(vf, M, is_oriented, reorient)

% plot_vf - plot a vector field with 
%   an optional image in the background.
%
% plot_vf(vf, M);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    M = [];
end
if nargin<3
    is_oriented = 1;
end
if nargin<4
    reorient = 0;
end

if reorient
    % reorient the vf to x>0
    epsi = sign(vf(:,:,1));
    I = find( epsi==0 );
    epsi(I) = 1;
    vf(:,:,1) = vf(:,:,1).*epsi;
    vf(:,:,2) = vf(:,:,2).*epsi;
end

n = size(vf,1);
p = size(vf,2);

x = 0:1/(n-1):1;
y = 0:1/(p-1):1;
[Y,X] = meshgrid(y,x);
 
hold on;
imagesc(x,y,M');
if is_oriented
    quiver(X,Y,vf(:,:,1),vf(:,:,2), 0.6);
else
    quiver(X,Y,vf(:,:,1),vf(:,:,2), 0.4);
    quiver(X,Y,-vf(:,:,1),-vf(:,:,2), 0.4);
end
% quiver(X,Y,-v(:,:,2),-v(:,:,1), 0.6);
axis xy;
axis equal;
axis off;
hold off;