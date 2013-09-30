function S = load_subimages_dataset(name,nbr,options)

% load_subimages_dataset - create a dataset of subimages.
%
%	S = load_subimages_dataset(name,nbr,options);
% or 
%   S = load_subimages_dataset(M,nbr,options);
%
%   S(:,:,i) is the ith sub-image of size options.dim of the full image M
%       (or loaded from given name).
%
%   M is an image.
%	name is a valid image name (eg. 'lena' or 'barb' if barb.gif or lena.png are in the path)
%	nbr is the number of images loaded
%	options.dim is the size of the extracted small squares
%	options.sampling is either 'rand' or 'uniform'.
%
%   Copyright (c) 2005 Gabriel Peyré

if ischar(name)
    M = load_image(name);
else
    M = name;
end

% force 1 channel (discard color)
M = sum(M,3);
n = size(M);

options.null = 0;

if isfield(options, 'dim')
    dim = options.dim;
else
    dim = [16 16];
end
S = zeros(dim(1),dim(2),nbr);

if isfield(options, 'sampling')
    sampling = options.sampling;
else
    sampling = 'rand';
end

if strcmp(sampling, 'rand')
    % random sampling
    x = 1 + floor( rand(nbr,1)*(n(1)-dim(1)) );
    y = 1 + floor( rand(nbr,1)*(n(2)-dim(2)) );
else
    % uniform sampling
    nbr = ceil( sqrt(nbr) )^2;
    x = linspace( 1, n(1)-dim(1)+1, sqrt(nbr) );
    y = linspace( 1, n(2)-dim(2)+1, sqrt(nbr) );
    [y,x] = meshgrid(y,x);
    x = x(:);
    y = y(:);
end

for i=1:nbr
    selx = x(i):x(i)+dim(1)-1;
    sely = y(i):y(i)+dim(2)-1;
    S(:,:,i) = M(selx,sely);
end


if isfield(options, 'equalize')
    equalize = options.equalize;
else
    equalize = 0;
end

if equalize
    % equalize mean and variance
    mu = sum(sum(S,1),2) / prod(dim);
    mu = repmat(mu,[dim(1) dim(2) 1]);
    S = S - mu;
    nu = sqrt( sum(sum(S.^2,1),2) / prod(dim) );
    nu = repmat(nu,[dim(1) dim(2) 1]);
    I = find(nu==0); nu(I) = 1;
    S = S ./ nu;
end