function plot_circle(q,r,options)

% plot_circle - display a collecion of circles
%
%   plot_circle(q,r,options);
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
draw_centers = getoptions(options, 'draw_centers', 1);
center_width = getoptions(options, 'center_width', 20);
npoints_circle = getoptions(options, 'npoints_circle', 60)+1;
color_circle = getoptions(options, 'color_circle', 'b');
color_center = getoptions(options, 'color_center', 'b');

if draw_centers
hh = plot(q(1,:), q(2,:), [color_center '.']);
end
set(hh, 'MarkerSize', center_width);
% draw circles
t = linspace(0,2*pi,npoints_circle);
for i=1:length(r)
    plot( sin(t)*r(i)+q(1,i), cos(t)*r(i)+q(2,i), color_circle );
end