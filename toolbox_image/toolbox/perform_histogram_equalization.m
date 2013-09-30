function x = perform_histogram_equalization(x,y,options)

% perform_histogram_equalization - perform histogram equalization
%
%   x = perform_histogram_equalization(x,y,options);
%
%   Change the values of x so that its ordered values match
%   the ordered values of y.
%
%   You can set
%       options.cols=1 to operate along columns
%       options.rows=1 to operate along rows
%       options.absval to operate only on absolute values.
%
%   You can also use
%       x = perform_histogram_equalization(x,'gaussian',options);
%       x = perform_histogram_equalization(x,'linear',options);
%   to impose gaussian or uniform density.
%   
%   Copyright (c) 2006 Gabriel Peyr?



options.null = 0;

if iscell(x)
    for i=1:length(x)
        x{i} = perform_histogram_equalization(x{i},y{i},options);
    end
    return;
end

if isfield(options ,'absval')
    absval = options.absval;
else
    absval = 0; 
end
if isfield(options ,'cols')
    cols = options.cols;
else
    cols = 0; 
end
if isfield(options ,'rows')
    rows = options.rows;
else
    rows = 0; 
end
if cols && rows
    error('You cannote specify both cols and rows');
end
if cols && size(x,2)>1
    if not(size(x,2)==size(y,2))
        error('x and y must have same number of columns');
    end
    for i=1:size(x,2)
        x(:,i) = perform_histogram_equalization(x(:,i),y(:,i),options);
    end
    return;
end
if rows && size(x,1)>1
    if not(size(x,1)==size(y,1))
        error('x and y must have same number of rows');
    end
    for i=1:size(x,1)
        x(i,:) = perform_histogram_equalization(x(i,:),y(i,:),options);
    end
    return;
end

sx = size(x);
x = x(:);

if isstr(y)
    switch y
        case 'gaussian'
            y = randn(length(x));
        case 'linear'
            y = linspace(0,1,length(x));
        otherwise 
            error('Unknown density.');
    end
end

y = y(:);


if absval
    s = sign(x);
    x = abs(x);
    y = abs(y);
end

[vx,Ix] = sort(x);
[vy,Iy] = sort(y);

nx = length(x);
ny = length(y);

ax = linspace(1,ny,nx);
ay = 1:ny;
vx = interp1(ay,vy,ax);
x(Ix) = vx;


if absval
    x = x .* s;
end

x = reshape(x,sx);