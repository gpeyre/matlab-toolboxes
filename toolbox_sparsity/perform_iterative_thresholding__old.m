function [X,E] = perform_iterative_thresholding(D,Y,options)

% perform_iterative_thresholding - perform L^1 minimization
%
%   [X,E] = perform_iterative_thresholding(D,Y,options);
%
%   Solve the L^1 penalized problem
%       min_X 1/2*|D*X-Y|_2^2 + lambda * |X|_1
%   It uses a proximal iteration
%           X <- Thresh_{mu*lambda}( X + mu*D'*( Y - D*X ) )
%   where mu is such that |D*z|^2<=mu^2*|z|^2
%
%   The number of iterations is options.niter.
%   options.thresh_type is 'hard' or 'soft'.
%
%   The threshold lambda=options.lambda control the degree of sparsity.
%   It can be modified linearly during the iterations by defining
%       options.lambda_min and options.lambda_max.
%   It can be estimated automatically using MAD estimator (useful for 
%       compressed sensing recovery).       
%
%   The dictionary D can be a 2D matrix, or it can be a callback function
%       D = @func;
%   where func has the form "y = func(D,x,dir,options)"
%   and dir=1 when computing y=D*x and dir=1 for y=D'*x.
%
%   You can provide an initial guess in options.X.
%
%   Copyright (c) 2006 Gabriel Peyre


options.null = 0;
niter           = getoptions(options, 'niter', 1);
thresh_type     = getoptions(options, 'thresh_type', 'soft');
verb            = getoptions(options, 'verb', 1);
strict_sparsity = getoptions(options, 'strict_sparsity', 10);
lambda          = getoptions(options, 'lambda', mean( sqrt( sum( Y.^2, 1 ) ) ) / 50);
lambda_min      = getoptions(options, 'lambda_min', lambda);
lambda_max      = getoptions(options, 'lambda_max', lambda_min*5);
mu_dampling         = getoptions(options, 'mu_dampling', 1);
tol = getoptions(options, 'tol', 0);
    
if strcmp(thresh_type, 'omp')
    if not(isfield(options, 'use_mex'))
        options.use_mex = 1;
    end
    options.nbr_max_atoms = strict_sparsity;
    options.verb = 0;
    X = perform_omp(D, Y, options);
    if issparse(X)
        X = full(X);
    end
    E = [];
    return;
end

e0 = norm(Y,'fro');

if niter>1
    options.niter = 1;
    E = [];
    if isnumeric(D) && not(isfield(options, 'mu'))
        [a,s,b] = svd(D); 
        options.mu = 1/max(diag(s))^2;
    end
    for i=1:niter
        options.lambda = lambda_max - (i-1)/(niter-1)*(lambda_max-lambda_min);
        if verb==1
            progressbar(i,niter);
        end
        X = perform_iterative_thresholding(D,Y,options);
        options.X = X;
        if tol>0 && norm( Y-D*X, 'fro')<tol*e0
            break;
        end
        if nargout>1
            if lambda_min>1e-4
                if isnumeric(D)
                    E(end+1) = 1/2*norm( Y-D*X, 'fro')^2 + lambda_min * sum( abs(X(:)) );
                else
                    E(end+1) = -1; % do not compute error in this case ...
                end
            else
                E(end+1) = sum( abs(X(:)) );
            end
        end
    end
    if nargout>1 && E(end)>E(1)
        warning('Iterative thresholding did not converge, you should lower options.mu_dampling.');
    end
    return;
end

if isfield(options, 'mu')
    mu = options.mu;
else
    if isnumeric(D)
        [a,s,b] = svd(D); 
        mu = 1/max(diag(s))^2;
    else
        warning('You should provide options.mu');
    end
end


if isfield(options, 'use_mad') && options.use_mad==1
    lambda_mca = 3;
    % estimate threshold via MAD estimation
    lambda = lambda_mca * mad(X,1)/0.6745;
end

if isfield(options, 'X') && not(isempty(options.X))
    X = options.X;
else
    if isnumeric(D)
        X = zeros(size(D,2), size(Y,2));
    else
        X = feval( D,  Y*0, -1, options );
    end
end

% sparse inversion
if isnumeric(D)
    X = X + mu * mu_dampling * D'*( Y - D*X );
else
    A = Y - feval( D, X, +1, options );
    X = X + mu * mu_dampling * feval( D, A, -1, options );
end

if strcmp(thresh_type, 'strict')
    t = options.strict_sparsity;
elseif strcmp(thresh_type, 'soft')
    t = lambda*mu_dampling*mu;
elseif strcmp(thresh_type, 'hard')
    % for hard thresholding, the scaling is different
    t = lambda*sqrt(mu_dampling*mu);
end
X = perform_thresholding( X, t, thresh_type);