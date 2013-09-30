function [x1,lun,err] = perform_douglas_rachford(A,y,options)

% perform_douglas_rachford - 
%
%   [x1,lun,err] = perform_douglas_rachford(A,y,options);
%
%   Solve noiseless Basis Pursuit
%       x1 = argmin_x \sum |x[k]|   subj.to A*x=y
%
%   A can be a matrix or an operator A(x,dir,options)
%   with dir=1 for A*x and dir=-2 for A^{+}*x (pseudo inverse).
%
%   options.x can be an initial guess
%
%   options.niter is the number of iterations
%
%   Special thanks to Jalal Fadili for useful helps on this algorithm.
%
%   Copyright (c) 2008 Gabirel Peyre

niter = getoptions(options, 'niter', 100);
mu = getoptions(options, 'mu', 3);

if isnumeric(A)
    pA = getoptions(options, 'pseudo_inv', []);
    if isempty(pA)
        pA = pinv(A);
    end
end

% initial guess
x = getoptions(options, 'x', [] );
if isempty(x)
    if isnumeric(A)
        x = pA*y;
    else
        x = feval(A, y, -2, options);
    end
end

% target result (error monitoring
x0 = getoptions(options, 'x0', [] );

drawiter = getoptions(options, 'drawiter',0);
verb = getoptions(options, 'verb', 1);
isreal = getoptions(options, 'isreal', 0);

solve_tv = getoptions(options, 'solve_tv', 0);

lun = [];
err = [];
if isnumeric(A)
    x1 = x + pA*(y-A*x);
else
    x1 = y - feval(A,x, +1, options);
    x1 = x + feval( A, x1, -2, options);
end
if drawiter
    clf;
end

nrefresh = 100;
ndisp = max(2, ceil(niter/nrefresh) );

for i=1:niter
    if verb
        progressbar(i,niter);
    end
    % compute the intermediary point
    %   x <- x-x1 + thresh( 2*x1-x )
    
    u = 2*x1-x;
    
%    sparse_weight = getoptions(options, 'sparse_weight', 0,1);
%    u = u./sparse_weight;
    if not(solve_tv)
        % L1 regularization
        u = perform_thresholding( u, mu, 'soft' );
    else
        % TV regularization
        opts.method = 'gradient';
        opts.lambda = mu;
        opts.niter = getoptions(options, 'niter_tv', 50);
        opts.niter_inner = 1;
        opts.display = 0;
        opts.verb = 0;
        [u1,err,tv,lalist,Err] = perform_tv_denoising(u,opts);
        % u1 = u - mu*div( xi, options )
%        opts.xi = grad( (u1-u)/mu, options );
        u = u1;
    end
%    u = u.*sparse_weight;
    
    x = x-x1 + u;
    
    % perform projection step on the constraints
    if isnumeric(A)
        x1 = x + pA*(y-A*x);
    else
        x1 = y - feval(A,x,+1, options);
        x1 = x + feval( A, x1, -2, options);
    end
    
    if isreal
 %       x1 = real(x1); % force real signal
    end
    % check for error decay
    if nargout>=3 && not(isempty(x0))
        err(end+1) = norm(x0-x1, 'fro')^2;
    end
    if nargout>=2
        if not(solve_tv)
            lun(end+1) = sum(abs(x1(:))); % mean(sum(abs(x1)));
        else
            lun(end+1) = compute_total_variation(x1,options);
        end
    end
    if drawiter && mod(i,ndisp)==1  % && ntest==1
        if nb_dims(x1)==1
            n = size(x1,1);
            plot( real(x1(:,1:min(end,3))) ); 
            axis([1,n,-1,1]); 
            axis tight;
            drawnow;
        else
            clf; imageplot(x1); drawnow;
        end
    end
end
% x1 = real(x1);
