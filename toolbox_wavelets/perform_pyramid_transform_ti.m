function MW = perform_pyramid_transform_ti( M, Jmin,options )

% perform_pyramid_transform_ti - gaussian pyramidal transform (1D or 2D).
%
% MW = perform_pyramid_transform_ti( M,Jmin,options )
%
%   'options' is a struct that can contains
%       - 'sigma': basis variance.
%       - 'bound': boundary extension (eiter 'sym' or 'per').
%
%   Decompose M as
%       M = MW{p} + sum_{k=1...p-1} MW{k}
%   with
%       MW{k} = M * ( G(sigma*2^j)-G(sigma*2^(j-1)) )
%   where G(s) is a gaussian filter with variance s.
%
%   Warning : this decomposition is redundant (it is
%   translation-invariant), and does not conserve energy at all ...
%   Good for computer vision applications ! 
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(M,1);
Jmax = log2(n)-1;

options.null = 1;
if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'sym';
end
if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 1;
end

if ~iscell(M)
    % FWD transform
    MW = {};
    M1 = M;
    for j = Jmax:-1:Jmin
        % compute the gaussian kernel
        s = sigma*2^(Jmax-j);   % width of the kernel in pixels
        p = min(s*4,n/2) + 4;  p = ceil(p/2)*2+1; % effective width
        g = compute_gaussian_filter( [p,p], [s/n,s/n], n);
        % filter
        M2 = perform_convolution(M,g,bound);
        MW{end+1} = M1-M2;
        M1 = M2;
    end
    MW{end+1} = M1;
else
    % BWD transform : just the sum (very redundant !)
    MW = M{end};
    for j=length(M)-1:-1:1
        MW = MW + M{j};
    end
end