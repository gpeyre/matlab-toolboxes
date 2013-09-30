function [H,X,Y] = compute_patch_library(M,w,options)

% [H,X,Y] = compute_patch_library(M,w,options);
%
%   M is the texture
%   w is the half-width of the patch, so that each patch
%       has size (2*w+1,2*w+1,s) where s is the number of colors.
%
%   H(:,:,:,i) is the ith patch (can be a color patch).
%   X(i) is the x location of H(:,:,:,i) in the image M.
%   Y(i) is the y location of H(:,:,:,i) in the image M.
%
%   options.sampling can be set to 'random' or 'uniform'
%   If options.sampling=='random' then you can define
%       options.nbr as the number of random patches and/or
%       option.locations_x and option.locations_y.
%
%   If w is a vector of integer size, then the algorithm 
%   compute a set of library, one for each size, 
%   and use the same location X,Y for each patch.
%   H{k} is the k-th library.
%
%   You can define options.wmax to avoid too big patch.
%   In this case, if w>wmax, the texture M will be filtered, and the patches
%   will be subsampled to match the size.
%
%   Copyright (c) 2006 Gabriel PeyrŽ

options.null = 0;
if length(w)>1
    % process the library of patches
    X = [];
    Y = [];
    for i=1:length(w)
        if ~isempty(X)
            option.locations_x = X;
            option.locations_y = Y;
            [H{i},X,Y] = compute_patch_library(M,w(i),options);
        end
    end
    return;
end

if isfield(options, 'sampling')
    sampling = options.sampling;
else
    sampling = 'random';
end
if isfield(options, 'wmax')
    wmax = options.wmax;
else
    wmax = 9;
end
wmax = min(w,wmax);

n = size(M,1);
n = size(M,1);
s = size(M,3);

if w>wmax
    % perform some smoothing to avoid aliasing when subsampling
    sigma = w/wmax; % in pixel size
    n1 = max( 16, round(n/8)*2+1 ); 
    h = create_gaussian_filter( [1 1]*n1,sigma/(2*n),[1 1]*n);
    M = perform_convolution(M,h);
end

% do some padding to avoid boundary problems
M = symmetric_extension(M,w);

ww = 2*w+1;
ww1 = 2*wmax+1;

if isfield(options, 'nbr') && strcmp(sampling, 'random')
    p = options.nbr;
else
    if strcmp(sampling, 'random')
        p = 2000;
    else
        p = n^2;
    end
end
    
if strcmp(sampling, 'random')
    if ~isfield(options, 'locations_x')
        X = w + 1 + floor(rand(p,1)*n);
        Y = w + 1 + floor(rand(p,1)*n);
    else
        X = options.locations_x(:)+w;
        Y = options.locations_y(:)+w;
        p = length(X);
    end
else
    [Y,X] = meshgrid(w+1:w+n, w+1:w+n);
    X = X(:);
    Y = Y(:);
end

H = zeros(ww1,ww1,s,p);
B = zeros(ww1,ww1,s);

if mod(w/wmax,1)==0    
    %%% in this case, a fast sampling can be used %%%
    [dY,dX] = meshgrid(-w:w/wmax:w,-w:w/wmax:w);
    Xp = repmat( reshape(X,[1,1,1,p]) ,[ww1 ww1 s 1]) + repmat(dX,[1 1 s p]);
    Yp = repmat( reshape(Y,[1,1,1,p]) ,[ww1 ww1 s 1]) + repmat(dY,[1 1 s p]);
    Cp = repmat( reshape(1:s,[1 1 s]), [ww1 ww1 1 p]);
    I = sub2ind(size(M), Xp,Yp,Cp);
    H = M(I);
    % remove padding in sampling location
    X = X-w; Y = Y-w;
    return;
end


%%% in this case, interpolation is needed %%%
xi = linspace(1,ww,ww1);
[Yi,Xi] = meshgrid(xi,xi);
h = waitbar(0,['Extracting patches of size ' num2str(w) ' ...']);
for i=1:p
    waitbar(i/p,h);
    x = X(i); y = Y(i);
    A = M( x-w:x+w, y-w:y+w,: );
    % need to perform interpolation
    for t=1:s
        B(:,:,t) = interp2(A(:,:,t),Xi,Yi)';
    end
    % perform sub-sampling
    H(:,:,:,i) = B;
end
close(h);

% remove padding in sampling location
X = X-w; Y = Y-w;