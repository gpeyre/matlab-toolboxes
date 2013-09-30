function plot_curve(c, M, str, line_width)

% plot_curve - plot a 2d curve.
% 
% plot_curve(c, M, str, line_width);
%
%   Note that 'c' can also be a cell array of curves.
%   - 'M' is a [n,p] image put in background (ranging in [0,1]x[0,1])
%   - 'c' is of size 2xm (ranging in [0,1]x[0,1])
%   - 'str' is the style of the curve.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    str = 'k';
end
if nargin<4
    line_width = 1;
end
if nargin<2
    M = [];
end

if iscell(c)
    if ~iscell(str)
        str1 = str; clear str;
        for i=1:length(c)
            str{i} = str1;
        end        
    end
    hold on;
    for i=1:length(c)
        plot_curve(c{i}, M, str{i}, line_width)
    end
    hold off;
    return;
end


if size(c,1)~=2
    c = c';
end
if size(c,1)~=2
    error('c should be of size 2xn');
end

hold on;

if ~isempty(M)
    [n,p] = size(M);
    x = 0:1/(n-1):1;
    y = 0:1/(n-1):1; % ( 0:1/(p-1):1 )* (p-1)/(n-1);
    % display image in X/Y frame
    M = M';
%    M = M(end:-1:1,:);
    imagesc(x,y,M);
    axis image;
end

plot( c(1,:), c(2,:), str, 'LineWidth', line_width );

axis tight;
hold off;