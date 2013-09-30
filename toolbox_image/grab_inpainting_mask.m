function [U,point_list] = grab_inpainting_mask(M, options)

% grab_inpainting_mask - create a mask from user input
%
%   U = grab_inpainting_mask(M, options);
%
%   Select set of point in an image (useful to select a region for
%   inpainting). The set of point is U==Inf.
%
%   options.r is the radius for selection (default r=5).
%
%   Selection stops with right click.
%
%   Set options.mode='points' to gather disconnected points.
%   Set options.mode='line' to gather connected lines.
%
%   Copyright (c) 2006 Gabriel Peyre


if nargin==3 && method==1
    U = grab_inpainting_mask_old(M, options);
    return;
end

options.null = 0;
r = getoptions(options, 'r', 5);
method = getoptions(options, 'mode', 'points');

if strcmp(method, 'line')
    if not(isfield(options, 'point_list'))
        [V,point_list] = pick_polygons(rescale(sum(M,3)),r);
    else
        point_list = options.point_list;
        V = draw_polygons(rescale(sum(M,3)),r,point_list);
    end        
    U = M; U(V==1) = Inf;
    return;
end

m = size(M,1);
n = size(M,2);
s = size(M,3);

U = sum(M,3)/3;
b = 1;
[Y,X] = meshgrid(1:n,1:m); 
point_list = [];
while b==1
    clf;
    hold on;
    imagesc(rescale(M)); 
    axis image; axis off; colormap gray(256);
    [y,x,b] = ginput(1);
    point_list(:,end+1) = [x;y];
    I = find((X-x).^2 + (Y-y).^2 <= r^2 );
    U(I) = Inf;
    for k=1:s
        Ma = M(:,:,k);
        Ma(I) = 0;
        M(:,:,k) = Ma;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sk = draw_polygons(mask,r,point_list)

sk = mask*0;
for i=1:length(point_list)
    pl = point_list{i};
    for k=2:length(pl)
        sk = draw_line(sk,pl(1,k-1),pl(2,k-1),pl(1,k),pl(2,k),r);
    end
end


function [sk,point_list] = pick_polygons(mask,r)

% pick_polygons - ask for the user to build a set of curves
%
%   sk = pick_polygons(mask,r);
%
%   mask is a background image (should be in [0,1] approx).
%
%   The user right-click on a set of point which create a curve.
%   Left click stop a curve.
%   Another left click stop the process.
%
%   Copyright (c) 2007 Gabriel Peyre


n = size(mask,1);

sk = zeros(n);
point_list = {};
b = 1;
while b(end)==1
    % draw a line
    clf;
    imagesc(mask+sk); axis image; axis off;
    colormap gray(256);
    [y1,x1,b] = ginput(1);
    pl = [x1;y1];
    while b==1
        clf;
        imagesc(mask+sk); axis image; axis off;
        [y2,x2,c] = ginput(1);
        if c~=1
            if length(pl)>1
                point_list{end+1} = pl;
            end
            break;
        end
        pl(:,end+1) = [x2;y2];
        sk = draw_line(sk,x1,y1,x2,y2,r);
        x1 = x2; y1 = y2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sk = draw_line(sk,x1,y1,x2,y2,r)


n = size(sk,1);
[Y,X] = meshgrid(1:n,1:n);
q = 80;
t = linspace(0,1,q);
x = x1*t+x2*(1-t); y = y1*t+y2*(1-t);
if r==0
    x = round( x ); y = round( y );
    sk( x+(y-1)*n ) = 1;
else
    for k=1:q
        I = find((X-x(k)).^2 + (Y-y(k)).^2 <= r^2 );
        sk(I) = 1;
    end
end





function U = grab_inpainting_mask_old(M, r)

% grab_inpainting_mask - create a mask from user input
%
%   U = grab_inpainting_mask(M, r);
%
%   r is the radius for selection (default r=5).
%
%   Selection stops with right click.
%
%   Copyright (c) 2006 Gabriel Peyr?

if nargin<2
    r = 5;
end


m = size(M,1);
n = size(M,2);
s = size(M,3);

U = sum(M,3)/3;
b = 1;
[Y,X] = meshgrid(1:n,1:m); 
while b==1
    clf;
    hold on;
    imagesc(rescale(M)); 
    axis image; axis off; colormap gray(256);
    [y,x,b] = ginput(1);
    I = find((X-x).^2 + (Y-y).^2 <= r^2 );
    U(I) = Inf;
    for k=1:s
        Ma = M(:,:,k);
        Ma(I) = 0;
        M(:,:,k) = Ma;
    end
end
