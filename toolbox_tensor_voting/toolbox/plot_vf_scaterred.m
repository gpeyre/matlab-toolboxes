function plot_vf_scaterred(vf, pos, is_oriented, reorient)

% plot_vf_scaterred - plot a vector field with 
%   an optional image in the background.
%
% plot_vf_scaterred(vf, pos, is_oriented, reorient);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    is_oriented = 1;
end
if nargin<4
    reorient = 0;
end


if reorient
    % reorient the vf to x>0
    epsi = sign(vf(2,:));
    I = find( epsi==0 );
    epsi(I) = 1;
    vf(1,:) = vf(1,:).*epsi;
    vf(2,:) = vf(2,:).*epsi;
end

% pack into matrix
n = size(vf,2);

p = ceil(sqrt(n));

U = zeros(p,p);
V = zeros(p,p);
X = zeros(p,p) + mean(pos(1,:));
Y = zeros(p,p) + mean(pos(2,:));
U(1:n) = vf(1,:);
V(1:n) = vf(2,:);
X(1:n) = pos(1,:);
Y(1:n) = pos(2,:);

 
hold on;
if is_oriented
    quiver(X,Y,U,V, 0.6);
else
    quiver(X,Y,U,V, 0.4);
    quiver(X,Y,-U,-V, 0.4);
end
% quiver(X,Y,-v(:,:,2),-v(:,:,1), 0.6);
axis xy;
axis equal;
axis off;
hold off;