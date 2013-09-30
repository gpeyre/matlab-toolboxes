function y = perform_windowed_dct_transform(x,w,q,n, options)

% perform_windowed_dct_transform - compute a local DCT transform
%
% Forward transform:
%   MF = perform_windowed_dct_transform(M,w,q,n, options);
% Backward transform:
%   M  = perform_windowed_dct_transform(MF,w,q,n, options);
%
%   w is the width of the window used to perform local computation.
%   q is the spacing betwen each window.
%
%   MF(:,:,,i,j) contains the transform around point ((i-1)*q,(j-1)*q)+1.
%
%   A typical use, for an redundancy of 2x2 could be w=2*q
%
%   Copyright (c) 2007 Gabriel Peyre

options.transform_type = 'dct';
y = perform_windowed_fourier_transform(x,w,q,n, options);


return;


%%%%%%%%%%% OLD CODE %%%%%%%%%%%%

options.null = 0;
if size(x,3)>1
    dir = -1;
    if nargin<4
        % assume power of 2 size
        n = q*size(x,3);
        n = 2^floor(log2(n));
    end
else
    dir = 1;
    n = size(x,1);
end

% perform sampling
t = 1:q:n-w+1;
[Y,X] = meshgrid(t,t);
p = size(X,1);

w = ceil(w/2)*2;
w1 = w/2;

t = 0:w-1;
[dY,dX] = meshgrid(t,t);

X = reshape(X,[1 1 p p]);
Y = reshape(Y,[1 1 p p]);
X = repmat( X, [w w 1 1] );
Y = repmat( Y, [w w 1 1] );
dX = repmat( dX, [1 1 p p] );
dY = repmat( dY, [1 1 p p] );

X1 = X+dX;
Y1 = Y+dY;

if 0
I = find(X1<1); X1(I) = 1-X1(I);
I = find(X1>n); X1(I) = 2*n+1-X1(I);
I = find(Y1<1); Y1(I) = 1-Y1(I);
I = find(Y1>n); Y1(I) = 2*n+1-Y1(I);
end


% build a weight function
if isfield(options, 'window_type')
    window_type = options.window_type;
else
    window_type = 'constant';
end

if strcmp(window_type, 'sin') || strcmp(window_type, 'sine')
    t = linspace(0,1,w);
    W = sin(t(:)*pi).^2;
    W = W * W';
elseif strcmp(window_type, 'constant')
    W = ones(w);
else
    error('Unkwnown winow.');
end

I = X1 + (Y1-1)*n;

if dir==1
    y = x(I) .* repmat( W, [1 1 p p] );
    y = my_dct_transform( y, +1 );
    
else
    x = my_dct_transform( x, -1 );
    weight = zeros(n); y = zeros(n);
    for i=1:p
    for j=1:p
        y(I(:,:,i,j)) = y(I(:,:,i,j)) + x(:,:,i,j);
        weight(I(:,:,i,j)) = weight(I(:,:,i,j)) + W;
    end
    end
    y = real( y./weight );
end


function y = my_dct_transform(x,dir)

for i=1:size(x,3)
    for j=1:size(x,4)
        y(:,:,i,j) = perform_dct_transform(x(:,:,i,j),dir);
    end
end

