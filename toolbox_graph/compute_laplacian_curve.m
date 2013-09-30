function [L,K] = compute_laplacian_curve(c,options)

% compute_laplacian_curve - compute a 2D laplacian
%
% [L,K] = compute_laplacian_curve(c,options)
%
%   c is a (n,d) curve of n points in R^d
%   options.is_closed==1 if the curve is closed.
%   options.method can be :
%       * 'combinatorial' : 0/1 weight
%       * 'exponential' : exp(-|xi-xj|^2/options.sigma^2)
%       * 'differential' : 1/(options.sigma + |xi-xj|)
%   options.normalization=1 to perform full normalization (ie the kernel K is symmetric)
%       otherwise, it will just sum to 1 (probabilistic matrix)
%
%   K is the weight matrix (symmetric positive)
%   L is the laplacian, L=(Id-K)/options.sigma.
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;

if isfield(options, 'is_closed')
    is_closed = options.is_closed;
else
    is_closed = 0;
end
if isfield(options, 'method')
    method = options.method;
else
    method = 'combinatorial';
end
if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.1;
end
if isfield(options, 'normalization')
    normalization = options.normalization;
else
    normalization = 1;
end


if size(c,1)<size(c,2)
    c = c';
end
n = size(c,1);

% compute the distance matrix
switch lower(method)
    case 'combinatorial'
        e = ones(n,1);
        D = spdiags([e e], [-1 1], n, n);
        if is_closed
            D(1,n) = 1;
            D(n,1) = 1;
        end
        I = find(D>0);
        D = D+Inf; D(I) = 0;

    otherwise
        % differential
        D = zeros(n) + Inf;
        for i=1:n
            x = c(i,:);
            if i<n
                a = i+1;
            else
                if ~is_closed
                    a = i-1;
                else
                    a = 1;
                end
            end
            if i>1
                b = i-1;
            else
                if ~is_closed
                    b = i+1;
                else
                    b = n;
                end
            end
            xa = c(a,:);
            xb = c(b,:);
            D(i,a) = sqrt( sum( (x-xa).^2 ) );
            D(i,b) = sqrt( sum( (x-xb).^2 ) );
        end
end

% compute the weight matrix
switch lower(method)
    case 'exponential'
        W = exp( -D.^2 / sigma^2 );
    case 'differential'
        W = 1. / (sigma+D);
    case 'combinatorial'
        W = double(D<sigma);
end

K = compute_diffusion_kernel(W,normalization);
K = sparse(K);
L = ( eye(n)-K )/sigma;