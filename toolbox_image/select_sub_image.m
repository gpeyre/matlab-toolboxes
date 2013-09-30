function [MM,sp,sr] = select_sub_image(M)

% select_sub_image - select a sub-part of an image.
%
%   [MM,sp,sr] = select_sub_image(M);
%
%   MM is the sub-image.
%   sp is the range of the image in pixels, ie. MM = M(sp(1):sp(2),sp(3):sp(4)).
%   sr is a scale of sr so that it fit into [0,1]².
%
%   Copyright (c) 2004 Gabriel Peyré

[n,p] = size(M);

imagesc(M);
axis image; axis off;
sp = getrect;

sp(1) = max(floor(sp(1)),1);    % xmin
sp(2) = max(floor(sp(2)),1);    % ymin
sp(3) = min(ceil(sp(1)+sp(3)),p);    % xmax
sp(4) = min(ceil(sp(2)+sp(4)),n);     % ymax

% swap X/Y
sp = [sp(2), sp(4), sp(1), sp(3)];
MM = M(sp(1):sp(2), sp(3):sp(4));

sr(1:2) = (sp(1:2)-1)/(n-1);
sr(3:4) = (sp(3:4)-1)/(p-1);