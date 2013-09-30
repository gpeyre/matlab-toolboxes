function y = callback_dictionary(x,dir,options)

% callback_atrou - MCA callback for matrix transform.
%
%   y = callback_dictionary(x,dir,options);
%
%   You must define options.D, options.w and options.n.
%
%   Copyright (c) 2007 Gabriel Peyre


if isfield(options, 'D')
    D = options.D;
else
    error('You must specify options.D');
end
if isfield(options, 'w')
    w = options.w;
else
    error('You must specify options.w');
end
if isfield(options, 'n')
    n = options.n;
else
    error('You must specify options.n');
end


% select the shift
if dir==1
    global shift_num;
    shift_num = floor( rand*w^2 ) + 1;
else
    global shift_num;
    if isempty(shift_num)
        error('shift undefined');
    end
end
[dY,dX] = meshgrid(0:w-1,0:w-1);
options.epsilon = [dX(shift_num) dY(shift_num)]; % offset

options.bound = 'per';
options.sub = w/2;

if isfield(options, 'sparse_coder')
    sparse_coder = options.sparse_coder;
else
    sparse_coder = 'omp';
    % sparse_coder = 'itthresh';
end

if dir==1
    s = size(x,3);
    % cut into patches
    H = compute_all_patch(x,w, options);
    p = size(H,4); % number of exemplar
    d = w^2*s; % dimensionality
    k = size(D,2); % number of atoms
    Y = reshape( H, d,p );
    % sparse coding
    switch sparse_coder
        case {'omp', 'mp'}
            X = perform_omp(D,Y,options);
            X = full(X);
        case 'itthresh'
            thresh_type = 'hard';
            % use iterative thresholding
            if isfield(options, 'T')
                lambda = options.T;
            else
                error('You must define options.T');
            end
            if isfield(options, 'niter_inversion')
                niter_inversion = options.niter_inversion;
            else
                niter_inversion = 6;
            end
            global Xprev;
            if isempty(Xprev)
                Xprev = zeros(k,p);
            end
            % conditioning of the dictionary
            [a,s,b] = svd(D); 
            mu = 1/max(diag(s))^2;
            % sparse inversion
            X = Xprev;
            for i=1:niter_inversion
                X = X + mu * D'*( Y - D*X );
                X = perform_thresholding( X, lambda*mu, thresh_type);
            end
            Xprev = X;
        otherwise
            error('Unknown coder.');
    end
    y = X;
else
    s = 1; % size(x,1)/w^2;
    % linear reconstruction
    Y = D*x;
    % patches
    H = reshape( Y, [w,w,s,size(Y,2)] );
    y = compute_all_patch(H, n, options);
end


if dir==-1
    shift_num = [];
end