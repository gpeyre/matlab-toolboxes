function plot_scattered(pos,z, resc)

% plot_scattered - plot 2D scattered data using triangulation.
%
%   plot_scattered(pos,z, resc);
%
%   if 'z' is ommited, it only plot a 2D scalar plot.
%   'resc==1' force a triangulation in [0,1]x[0,1].
%
%   'pos' should be of size 3xn or 2xn
%
%   Copyright (c) 2004 Gabriel Peyré

% eventually swap
if size(pos,1)~=2 && size(pos,1)~=3
    pos = pos';
end
if size(pos,1)~=2 && size(pos,1)~=3
    error('pos is not of correct format.');
end
if nargin<2
    z = [];
end

if size(pos, 1)==3
    % plot3( pos(1,:), pos(2,:), pos(3,:), '.' );
    if ~isempty(z)
        scatter3(pos(1,:),pos(2,:),pos(3,:),12,z,'filled');
    else
        scatter3(pos(1,:),pos(2,:),pos(3,:),12,pos(1,:),'filled');
    end

    axis tight;
    axis equal;
    cameramenu;

    return;
end

if nargin<3
    if ~isempty(z)
        scatter(pos(1,:),pos(2,:),12,z,'filled');
    else
        scatter(pos(1,:),pos(2,:),12,'filled');
    end
    % plot( pos(1,:), pos(2,:), '.' );    
%    axis equal;
    axis tight;
    return;
end

% 2D + height plot

if nargin<3
    resc = 0;
end

X = pos(1,:)';
Y = pos(2,:)';

if resc>0
    X = rescale(X)*resc;
    Y = rescale(Y);
end

TRI = delaunay(X,Y);

trisurf(TRI,X,Y,z);
view(45,20);

axis tight;
axis equal;
cameramenu;