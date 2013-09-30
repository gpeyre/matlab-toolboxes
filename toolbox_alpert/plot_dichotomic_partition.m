function plot_dichotomic_partition(X, part, options)

% plot_dichotomic_partition - plot clustered points in R^2 or R^3.
%
%   plot_dichotomic_partition(X, part);
%
%   Copyright (c) 2004 Gabriel Peyré

if size(X,1)>size(X,2)
    X = X';
end

options.null = 0;

if isfield(options, 'point_size')
    point_size = options.point_size;
else
    point_size = 25;
end


d = size(X,1);
if d>3
    d = 3;
    X = X(1:3,:);
end

hold on;
for i=1:length(part)
    I = part{i};
    if d==2
        scatter(X(1,I),X(2,I), point_size, 'filled');
    elseif d==3
        scatter3(X(1,I),X(2,I),X(3,I), point_size, 'filled');
    end
end
hold off;