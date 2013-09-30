function M = play_movie(A, nbr, fps)

% play_movie - play a volumetric data as an animation
%
%   play_movie(A);  % default
%   M = play_movie(A, nbr, fps);
%
%   A is an (n,n,nbr_frame) matrix.
%
%   Copyright (c) 2006 Gabriel Peyré


display_nbr = 1;

if nargin<2
    nbr = Inf;
end
if nargin<3
    fps = 24;
end

p = size(A,3);

A = rescale(A);
n = size(A,1);


fig=figure;
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
    'NextPlot','replace','Visible','off')
mov = avifile('example.avi');

x = -pi:.1:pi;
radius = [0:length(x)];
for i=1:length(x)
    h = patch(sin(x)*radius(i),cos(x)*radius(i),[abs(cos(x(i))) 0 0]);
    set(h,'EraseMode','xor');
    F = getframe(gca);
    mov = addframe(mov,F);
end
mov = close(mov);

close all;
for i=1:p
    clf;
    image(A(:,:,i)*256);
    colormap gray(256);
    axis image; axis off;
    if display_nbr
        text(1,1,num2str(i));
    end
    M(i) = getframe;
    % M(i) = im2frame( A(:,:,i), map );
end

if nargout==0
    movie(M,fps);
end