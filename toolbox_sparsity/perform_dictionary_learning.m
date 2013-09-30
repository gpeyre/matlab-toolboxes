function [D,X,E] = perform_dictionary_learning(Y,options)

% perform_dictionary_learning - learn a dictionnary using K-SVD algorithm
%
%   [D,X,E] = perform_dictionary_learning(Y,options)
%
%   Y is a matrix of size (n,m) containing m column examplar
%       vector in R^n.
%
%   D is a dictionnary matrix of size (n,K) of K vectors in R^n
%       that should approximate well the vectors of Y 
%       with few components.
%       K is given by options.K
%
%   The algorihtm perform the folowing optimization jointly on D and X
%       min_{D,X}  |Y-D*X|^2   subject to |X_k|_0<=s
%           subject to columns of D having unit norm.
%   where s=options.nbr_max_atoms.
%
%   It does this optimisation using block coordinate descent.
%       * Optimization of X knowing D amount to sparse coding, solved using
%       iterative thresholding:
%           X <- Thresh_lambda( X + D'*( Y - D*X ) )
%       This step is repeated options.niter_inversion times.
%       * Optimization of D knowing X amount to L2 best fit:
%           D <- Y*X^+   where   X^+ = X'*(X*X')^{-1}
%
%   This algorithm is very much inspired by the MOD algorithm
%   for dictionary design.
%
%   The number of iteration of the algorithm is given in
%       options.niter_learning.
%   The number of iteration for sparse coding during each step of learning
%       is given by options.niter_inversion.
%   The sparse coder used is set in options.sparse_coder to either 'omp' or
%   'mp' or 'itthresh'.
%
%   options.learning_method can be set to:
%       'ksvd': in this case, the dictionary update is computed using the
%           K-SVD algorithm explained in 
%           "K-SVD: An Algorithm for Designing Overcomplete Dictionaries
%                   for Sparse Representation"
%           Michal Aharon Michael Elad Alfred Bruckstein, 2006
%       'mod': in this case, the dictionary update is computed using the L2
%           best fit as proposed by
%           "Method of optimal directions for frame design",
%           K. Engan, S. O. Aase, and J. H. Hus?y,
%           in Proc. ICASSP, Phoenix, AZ, Mar. 1999, pp. 2443?2446. 
%       'randomized': apply randomly one or the other algorithm.
%       'modortho' same as mod, but the dictionary is constrained to be
%           orthogonal. In this case, the dictionary must be square. The
%           idea is taken from           
%               "Learning Unions of Orthonormal Bases with Thresholded Singular Value Decomposition"
%               S Lesage, R Gribonval, F Bimbot, and L Benaroya
%               in Proc. ICASSP'05, vol. V, p. 293--296
%       'grad': same as MOD but with a pojected gradient descent.
%
%   If you work with missing data, then you need to provide
%   mask=options.mask which should a matrix of the same size as Y and
%   mask(i,j)=0 whenever an entry Y(i,j) is missing (should be Y(i,j)=0)
%   and mask(i,j)=1 for valid data.
%
%   MOD is usualy faster but KSVD gives a better optimization of the energy.
%   KSVD is very efficient when a very low sparsity is required (set options.nbr_max_atoms to a small constant).
%
%   Copyright (c) 2007 Gabriel Peyre


[n,m] = size(Y);

options.null = 0;
niter           = getoptions(options, 'niter_learning', 10);
K               = getoptions(options, 'K', n);
niter_inversion = getoptions(options, 'niter_inversion', 4);
lambda          = getoptions(options, 'lambda', mean( sqrt( sum( Y.^2, 1 ) ) ) / 50);
sparse_coder    = getoptions(options, 'sparse_coder', 'omp');
mu              = getoptions(options, 'mu', 1/K);
init_dico       = getoptions(options, 'init_dico', 'input');
mu_dampling     = getoptions(options, 'mu_dampling', 1);
lambda_min      = getoptions(options, 'lambda_min', lambda);
lambda_max      = getoptions(options, 'lambda_max', lambda_min*4);
mask            = getoptions(options, 'mask', []);
use_bootstrapping = getoptions(options, 'use_bootstrapping', 0);
verb            = getoptions(options, 'verb', 1);
learning_method = getoptions(options, 'learning_method', 'mod');

if not(isfield(options, 'nbr_max_atoms'))
    warning('You should set the sparsity options.nbr_max_atoms.');
    options.nbr_max_atoms = 4;
end

if niter>1
    options.niter_learning = 1;
    E = [];
    for i=1:niter
        options.lambda = lambda_max - (i-1)/(niter-1)*(lambda_max-lambda_min);
        if verb
            progressbar(i,niter);
        end
        [D,X] = perform_dictionary_learning(Y,options);
        options.D = D;
        options.X = X;
        if strcmp(sparse_coder, 'strict') || strcmp(sparse_coder, 'omp') || strcmp(sparse_coder, 'mp')
            if isempty(mask)
                E(end+1) = norm( Y-D*X, 'fro');
            else
                E(end+1) = norm( (Y-D*X).*mask, 'fro');                
            end
        else
            E(end+1) = 1/2*norm( Y-D*X, 'fro')^2 + lambda_min * sum( abs(X(:)) );
        end
    end
    % sort the basis vector according to energy
    e = sum( X.^2, 2 );
    [tmp,I] = sort(e); I = I(end:-1:1);
    D = D(:,I); X = X(I,:);
    return;
end

if isfield(options, 'D') && not(isempty(options.D))
    D = options.D;
else
    % intialize dictionary
    switch init_dico
        case 'input'
            sel = randperm(m); sel = sel(1:K);
            D = Y(:,sel);
        case 'haar1d'
            D = compute_haar_matrix(n,1);
        case 'haar2d'
            D = compute_haar_matrix([sqrt(n) sqrt(n)],1);
        case 'rand'
            D = randn(n,K);
    end
    if size(D,2)<K
        D = [D; randn(n,K-size(D,2))];
    elseif size(D,2)>K
        D = D(:,1:K);
    end
    D = D - repmat( mean(D), n,1 );
    D(:,1) = 1; % enforce low pass
    D = D ./ repmat( sqrt(sum(D.^2,1)), n,1 );
end

if isfield(options, 'X') && not(isempty(options.X))
    X = options.X;
else
    X = zeros(size(D,2), size(Y,2));
end


% sparse inversion
if isempty(mask)
    if not(strcmp(sparse_coder, 'omp')) && not(strcmp(sparse_coder, 'mp'))
        % sparse code wiht iterative thresholdings
        % compute condition number
        if isfield(options, 'mu')
            mu = options.mu;
        else
            [a,s,b] = svd(D);
            mu = 1/max(diag(s))^2;
        end
        E = [];
        for i=1:niter_inversion
            X = X + mu * mu_dampling * D'*( Y - D*X );
            if strcmp(sparse_coder, 'soft')
                t = lambda*mu_dampling*mu;
            elseif strcmp(sparse_coder, 'hard')
                % for hard thresholding, the scaling is different
                t = lambda*sqrt(mu_dampling*mu);
            end
            X = perform_thresholding( X, t, sparse_coder);
            E(end+1) = 1/2*norm( Y-D*X, 'fro')^2 + lambda * sum( abs(X(:)) );
        end
        if E(end)>E(1)
            warning('Iterative thresholding did not converge, you should lower options.mu_dampling.');
        end
    else
        % sparse code with matching pursuit
        if not(isfield(options, 'use_mex'))
            options.use_mex = 0;
        end
        if strcmp(learning_method, 'modortho')
            M = options.nbr_max_atoms; % number of coefficients to keep
            p = size(Y,2); % number of signals
            d = size(Y,1);
            % use simple best k-term in ortho-basis
            X = D'*Y; % coefficients
            [tmp,I] = sort(abs(X));
            I = I(end-M,:) + (0:p-1)*d;
            T = repmat( abs(X(I)), [d 1] );
            X = X.*(abs(X)>T);
        else
            options.verb = 0;
            X = perform_omp(D, Y, options);
            if issparse(X)
                X = full(X);
            end
        end
    end
else
    if use_bootstrapping
        % fill-in the missing values
        Y0 = Y; Y = D*X; 
        Y(mask==1) = Y0(mask==1);
        % sparse code
        % options.learning_method = 'mod';
        options.verb = 0;
        X = full( perform_omp(D, Y, options) );
    else
        % sparse code with atom with missing information
        for k=1:size(Y,2)
            % progressbar(k,size(Y,2));
            I = find(mask(:,k));
            X(:,k) = full( perform_omp( D(I,:), Y(I,k), options) );
        end
        options.learning_method = 'grad';
    end
end


if strcmp(learning_method, 'randomized')
    randomized_ksvd_proportion = getoptions(options, 'randomized_ksvd_proportion', .3);    
    if rand<randomized_ksvd_proportion
        learning_method = 'ksvd';
    else
        learning_method = 'mod';
    end
end

if strcmp(learning_method, 'ksvd')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       KSVD algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = 1e-3;
    % the first atoms is supposed to be constant
    for k=2:K
        % find the exemplar that are using dictionnary basis k
        I = find( abs(X(k,:))>T );
        if ~isempty(I)
            % compute redisual
            if 0
                D0 = D; D0(:,k) = 0;
                E0 = Y - D0*X;
                % restrict to element actually using k
                E0 = E0(:,I);
            else
                S = X(:,I);
                S(k,:) = 0;
                E = Y(:,I) - D*S;
            end
            % perform SVD
            [U,S,V] = svd(E);
            D(:,k) = U(:,1);
            X(k,I) = S(1) * V(:,1)';
        end
    end
elseif strcmp(learning_method, 'mod')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       MOD algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dictionary fit
    warning off;
    % D = Y * X'*(X*X')^(-1);
    D = Y * pinv(X);
    warning on;
elseif strcmp(learning_method, 'grad')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Gradient descent of the MOD energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(mask)
        mask = Y*0+1;
    end
    if isfield(options, 'niter_grad')
        niter_grad = options.niter_grad;
    else
        niter_grad = 30;
    end
    if isfield(options, 'lambda_grad')
        lambda_grad = options.lambda_grad;
    else
        lambda_grad = .01;
    end
    % dictionary update by gradient descent
    Res = (Y-D*X).*mask;
    errg = [];
    for it=1:niter_grad
        D = D + lambda_grad*Res*X';
        % normalize
        D = normalize_dictionary(D);
        % error
        Res = (Y-D*X).*mask;
        errg(end+1) = norm(Res,'fro');
    end
    if errg(end)>errg(1)
        warning('Gradient descent did not converged');
    end
elseif strcmp(learning_method, 'modortho')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Orthogonal MOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n~=K
        error('For modortho learning, n=K is required.');
    end
    [U,S,V] = svd(Y*X');
    D = U*V';
else
    error('Unknown learning method');
end
D = normalize_dictionary(D);
E = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = normalize_dictionary(D)
n = size(D,1);
D = D - repmat( mean(D),[n 1] );
D(:,1) = 1; % enforce low pass
% Normalize
d = sqrt(sum(D.^2,1));
I = find(d<1e-9); d(I) = 1;
D(:,I) = randn(size(D,1),length(I));
d = sqrt(sum(D.^2,1));
D = D ./ repmat( d, n,1 );