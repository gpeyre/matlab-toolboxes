function display_image_layout(Mlist, TitleList, a, b)

% display_image_layout - display a set of images together
% 
%   display_image_layout(Mlist, TitleList, a, b);
%
%   Image should be in [0,1].
%   b is the number of row (optional)
%   a is the number of columns (optional)
%
% Example:
%   display_image_layout( {A B}, {'image 1' 'image 2'} );
%
%   Copyright (c) 2007 Gabriel Peyre

m = length(Mlist);

if nargin<3
    a = round( sqrt(m)/1.2 );
end
if nargin<4
    b = ceil( m/a );
end
if a*b<m
    warning('You must increase a and b');
    m = a*n;
end

clf;
ax = [];
for i=1:m
    M = Mlist{i};
    ax(end+1) = subplot(a,b,i);
    image( clamp(M)*255 ); axis image; axis off;
    title(TitleList{i});
end

linkaxes(ax, 'xy');
colormap gray(256);