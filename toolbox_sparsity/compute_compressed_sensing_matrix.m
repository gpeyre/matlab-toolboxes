function M = compute_compressed_sensing_matrix(p,n, type, options)

% compute_compressed_sensing_matrix - compute a CS matrix
%
%   M = compute_compressed_sensing_matrix(p,n, type, options);
%
%   M is a (p,n) matrix.
%   type can be
%       'gaussian': gaussian entries
%       'randnormed': random unit norm columns
%       'bernouilli': random +-1/sqrt(p)
%       'orthoproj': random projector M'*M=Id
%       'fourier': random Fourier rows
%       'sincos': random cosine/sine rows
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<3
    type = 'gaussian';
end

options.null = 0;
switch lower(type)
    case 'gaussian'
    	M = randn(p,n) / sqrt(p);
    case {'randnormed' 'rand'}
        M = randn(p,n);
        M = M ./ repmat( sqrt(sum(M.^2)), [p 1] );
    case 'bernouilli'
        M = sign(randn(p,n)) / sqrt(p);
    case {'orthoproj' 'randproj'}
        M = randn(p,n);
        [M,R] = qr(M'); M = M(:,1:p)'; 
    case 'fourier'
        [X,Y] = meshgrid(0:n-1:0:n-1);
        M = exp( 2i*pi/n * X.*Y ) / sqrt(p);
        sel = randperm(n); sel = sel(1:p);
        M = M(sel,:);
    case 'sincos'
        [X,Y] = meshgrid(0:n-1,0:n/2);
        M = cos( 2*pi/n * X.*Y );
        [X,Y] = meshgrid(0:n-1,1:n/2-1);
        M = [M; sin( 2*pi/n * X.*Y )];
        M = M ./ repmat( sqrt(sum(M.^2,2)), [1 n] );
        sel = randperm(n); sel = sel(1:p);
        M = M(sel,:) * sqrt( n/p );
    case 'projoptim'
        niter = getoptions(options, 'niter_projoptim', 40);
        M = randn(p,n);
        for i=1:niter
            % unit normed
            M = M ./ repmat( sqrt(sum(M.^2)), [p 1] );
            % projector
            [U,S,V] = svd(M, 'econ');
            S = diag(diag(S)*0+1);
            M = U*S*V';            
        end
    otherwise 
        error('Unknown matrice type');
end

return;
%% old code

if nargin<4
    normtype = 'normalize';
end

switch lower(type)
    case 'random'
        M = randn(n,m);
    case 'bumps' 
        normtype = 'none';
        % m bumps
        sigma = 0.05;
        x = linspace(0,1,m+1)'; x(end) = []; % bump centers
        M = zeros(n,m);
        t = linspace(0,1,n+1)'; t(end) = [];
        for i=1:m
            u = t-x(i); u(u>1/2) = 1-u(u>1/2);
            M(:,i) = exp( -u.^2 / (2*sigma^2)  );
        end
end

switch lower(normtype)
    case 'normalize'
        % lines should be of unit norm
        d = sqrt( sum(M.^2,1) );
        M = M ./ repmat( d, [size(M,1) 1] );
    case 'orthogonalize'
        % M*M' = Id
        M = orth(M')';
    case 'none'
    otherwise
        error('Unknown normalization');
end