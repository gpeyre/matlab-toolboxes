function [M,eOrth,eAngle,eNorm] = perform_cs_matrix_optimization(M, options)

% perform_cs_matrix_optimization - enhance sensing matrix
%
%   Enforces orthogonality of rows and unit norm of columns of the sensing
%   matrix M.
%
%   This tends to give better CS results.
%
%   Copyright (c) 2007 Gabriel Peyre


n = size(M,1);
m = size(M,2);

if isfield(options, 'use_norm')
    use_norm = options.use_norm;
else
    use_norm = 1;
end
if isfield(options, 'use_angle')
    use_angle = options.use_angle;
else
    use_angle = 0;
end
if isfield(options, 'use_orth')
    use_orth = options.use_orth;
else
    use_orth = 1;
end
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 20;
end

eOrth = [];
eAngle = [];
eNorm = [];
lambda = 0.05;
for i=1:niter
    progressbar(i,niter);
    if use_orth
        % orthogonality constraint
        eOrth(end+1) = norm( M*M' - eye(n)*m/n );
        % M = orth(M')';
        [U,S,V] = svd(M);
        S(S>0) = 1;
        M = U*S*V';
        M = M(1:n,:) * sqrt(m/n);
    end
    % norm constraint
    if use_norm
        d = sqrt( sum( M.^2, 1 ) );
        eNorm(end+1) = sum( abs(d-1) );
        M = M./repmat( d, [n 1] );
    end
    if use_angle
        % angle constraint gradient descent
        D = M'*M-eye(m);
        eAngle(end+1) = sum( sum(D.^2) );
        nitergrad = 10;
        for i=1:nitergrad
            D = M'*M-eye(m);
            G = 2 * M * D; %  ( D + D' );
            M = M - lambda * G;
        end
    end
end