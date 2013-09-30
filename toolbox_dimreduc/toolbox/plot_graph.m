function plot_graph(A,xy, col)

% plot_graph - display a 2D or 3D graph.
%
% plot_graph(A,xy, col);
%
%   Copyright (c) 2006 Gabriel PeyrŽ

if size(xy,1)>size(xy,2)
    xy = xy';
end

if nargin<3
    col = 'k.-';
end

if size(xy,1)==2
    if ~isstr(col)
        col = [];
    end
    % 2D display
    gplot(A,xy', col);
elseif size(xy,1)==3
    % 3D display
    if isstr(col)
        gplot3(A,xy', col);
    else
        hold on;
        gplot3(A,xy', 'k');
        plot_scattered(xy, col);
        hold off;
        view(3);
    end
else
    error('Works only for 2D and 3D graphs');
end


function gplot3(A,xy,lc)


[i,j] = find(A);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-separated list of line segments,
% rather than individual segments.

X = [ xy(i,1) xy(j,1) repmat(NaN,size(i))]';
Y = [ xy(i,2) xy(j,2) repmat(NaN,size(i))]';
Z = [ xy(i,3) xy(j,3) repmat(NaN,size(i))]';
X = X(:);
Y = Y(:);
Z = Z(:);

if nargin<3,
    plot3(X, Y, Z)
else
    plot3(X, Y, Z, lc);
end