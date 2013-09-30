function y = perform_windowed_fourier_transform(x,w,q,n, options)

% perform_windowed_fourier_transform - compute a local Fourier transform
%
% Forward transform:
%   MF = perform_windowed_fourier_transform(M,w,q,n, options);
% Backward transform:
%   M  = perform_windowed_fourier_transform(MF,w,q,n, options);
%
%   w is the width of the window used to perform local computation.
%   q is the spacing betwen each window.
%
%   MF(:,:,,i,j) contains the spectrum around point ((i-1)*q,(j-1)*q)+1.
%
%   A typical use, for an redundancy of 2x2 could be w=2*q+1
%
%   options.bound can be either 'per' or 'sym'
%
%   options.normalization can be set to
%       'tightframe': tight frame transform, with energy conservation.
%       'unit': unit norm basis vectors, usefull to do thresholding
%
%   Copyright (c) 2006 Gabriel Peyre

force_energy = 1;

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

bound = getoptions(options, 'bound', 'sym');
transform_type = getoptions(options, 'transform_type', 'fourier');
normalization = getoptions(options, 'normalization', 'tightframe');

% perform sampling
t = 1:q:n+1;
[Y,X] = meshgrid(t,t);
p = size(X,1);

if mod(w,2)==1
% w = ceil((w-1)/2)*2+1;
    w1 = (w-1)/2;
    t = -w1:w1;
else
    t = -w/2+1:w/2;
end
[dY,dX] = meshgrid(t,t);

X = reshape(X,[1 1 p p]);
Y = reshape(Y,[1 1 p p]);
X = repmat( X, [w w 1 1] );
Y = repmat( Y, [w w 1 1] );
dX = repmat( dX, [1 1 p p] );
dY = repmat( dY, [1 1 p p] );

X1 = X+dX;
Y1 = Y+dY;

switch lower(bound)
    case 'sym'
        X1(X1<1) = 1-X1(X1<1);
        X1(X1>n) = 2*n+1-X1(X1>n);
        Y1(Y1<1) = 1-Y1(Y1<1);
        Y1(Y1>n) = 2*n+1-Y1(Y1>n);
    case 'per'
        X1 = mod(X1-1,n)+1;
        Y1 = mod(Y1-1,n)+1;
end
        

% build a weight function
window_type = getoptions(options, 'window_type', 'sin');
eta = getoptions(options, 'eta', 1);

if strcmp(window_type, 'sin')
    t = linspace(0,1,w);
    W = sin(t(:)*pi).^2;
    W = W * W';
elseif strcmp(window_type, 'constant')
    W = ones(w);
else
    error('Unkwnown winow.');
end

I = X1 + (Y1-1)*n;

%% renormalize the windows
weight = zeros(n);
for i=1:p
    for j=1:p
        if force_energy
            % special care ... because of overlaping
            for s=1:size(I,1)
                for t=1:size(I,2)
                    weight(I(s,t,i,j)) = weight(I(s,t,i,j)) + W(s,t).^2;
                end
            end
        else
            weight(I(:,:,i,j)) = weight(I(:,:,i,j)) + W.^2;
        end
    end
end
weight = sqrt(weight);
Weight = repmat(W, [1 1 p p]);
for i=1:p
    for j=1:p
        if force_energy
            % special care ... because of overlaping
            for s=1:size(I,1)
                for t=1:size(I,2)
                    Weight(s,t,i,j) = Weight(s,t,i,j) ./ weight(I(s,t,i,j));
                end
            end
        else
            Weight(:,:,i,j) = Weight(:,:,i,j) ./ weight(I(:,:,i,j));            
        end
    end
end


if strcmp(normalization, 'unit')
    if strcmp(transform_type, 'fourier')
        % for Fourier it is easy
        Renorm = sqrt( sum( sum( Weight.^2, 1 ), 2 ) )/w;
    else
        % for DCT it is less easy ...
        % take a typical window in the middle of the image
        weight = Weight(:,:,round(end/2),round(end/2));
        % compute diracs
        [X,Y,fX,fY] = ndgrid(0:w-1,0:w-1,0:w-1,0:w-1);
        A = 2 * cos( pi/w * ( X+1/2 ).*fX ) .* cos( pi/w * ( Y+1/2 ).*fY ) / w;
        A(:,:,1,:) = A(:,:,1,:) / sqrt(2); % scale zero frequency
        A(:,:,:,1) = A(:,:,:,1) / sqrt(2); 
        A = A .* repmat( weight, [1 1 w w] );
        Renorm = sqrt( sum( sum( abs(A).^2, 1 ),2  ) );
    end
    Renorm = reshape( Renorm, w,w );
end
    
%% compute the transform
if dir==1
    y = zeros(eta*w,eta*w,p,p);
    if mod(w,2)==1
        m = (eta*w+1)/2; w1 = (w-1)/2;
        sel = m-w1:m+w1;
    else
        m = (eta*w)/2+1; w1 = w/2;
        sel = m-w1:m+w1-1;
    end
    y(sel,sel,:,:) = x(I) .* Weight;
    % perform the transform
    y = my_transform( y, +1, transform_type );
    % renormalize if necessary
    if strcmp(normalization, 'unit')
        y = y ./ repmat( Renorm, [1 1 p p] );
    end
else
    if strcmp(normalization, 'unit')
        x = x .* repmat( Renorm, [1 1 p p] );
    end
    x = my_transform( x, -1, transform_type );
    x = real( x.*Weight );
    y = zeros(n);
    for i=1:p
        for j=1:p            
            if force_energy
                % special care ... because of overlaping
                for s=1:size(I,1)
                    for t=1:size(I,2)
                        y(I(s,t,i,j)) = y(I(s,t,i,j)) + x(s,t,i,j);
                    end
                end
            else
                y(I(:,:,i,j)) = y(I(:,:,i,j)) + x(:,:,i,j);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_transform(x,dir,transform_type)

% my_transform - perform either FFT or DCT with energy conservation.
%   Works on array of size (w,w,a,b) on the 2 first dimensions.

w = size(x,1);
if strcmp(transform_type, 'fourier')
    % normalize for energy conservation
    if dir==1
        y = myfftshift( fft2(x) ) / w;
    else
        y = ifft2( myifftshift(x*w) );
    end
elseif strcmp(transform_type, 'dct')
    for i=1:size(x,3)
        for j=1:size(x,4)
            y(:,:,i,j) = perform_dct_transform(x(:,:,i,j),dir);
        end
    end
else
    error('Unknown transform');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = myfftshift(x)
x = fftshift(x,1);
x = fftshift(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = myifftshift(x)
x = ifftshift(x,2);
x = ifftshift(x,1);

