function [x,err,lun,Tlist,Err] = perform_iterative_thresholding(U,y,options)

% perform_iterative_thresholding - perform L^1 minimization
%
%   [x,err,lun,Tlist] = perform_iterative_thresholding(U,y,options);
%
%   Solve the L^1 penalized problem
%       min_x 1/2*|U*D*x-y|_2^2 + T * |X|_1
%   It uses a proximal iteration
%           x <- Thresh_{T/tau}( x + 1/tau * (U*D)'*( y - (U*D)*x ) )
%   where tau is such that |D*U*z|^2 <= tau*|z|^2
%
%   If T==0, then solve
%       min_x |x|_1   subject to U*D*x=y
%
%   Usually, U is the operator to invert, and D is a sparsifying
%   dictionary. D is optional (it is given in options.D).
%
%   The number of iterations is options.niter.
%   options.thresh is 'hard' or 'soft'.
%
%   The threshold T=options.T control the degree of sparsity.
%   It can be modified linearly during the iterations by defining
%       options.Tmin and options.Tmax. 
%
%   In order to solve for
%       min_x |X|_1   subject to    |U*D*x-y|_2 <= etgt
%   you can set the value of options.etgt. In this case, the algorithm will
%   update the value of T through the iteration to match the residual error
%   provided.
%
%   The dictionary D can be a 2D matrix, or it can be a callback function
%       D = @func;
%   where func has the form "y = func(x,dir,options)"
%   and dir=+1 when computing y=D*x (reconstruction), 
%   and dir=-1 when computing x=D'*y, and dir=-2 for x=D^+*y
%   (pseudo inverse).
%
%   Same thing for U.
%
%   You can provide an initial guess in options.x.
%   You can provide an oracle solution in options.x0.
%
%   Copyright (c) 2006 Gabriel Peyre


niter = getoptions(options, 'niter', 200);
D = getoptions(options, 'D', []);

thresh = getoptions(options, 'thresh','soft');

T = getoptions(options, 'T', []);
Tmax = .1;
Tmin = 0;
Tmax = getoptions(options, 'Tmax', Tmax);
Tmin = getoptions(options, 'Tmin', Tmin);
if not(isempty(T))
    Tmax = T; Tmin = T;
end
Tlist = linspace(Tmax,Tmin,niter);

exlude_low = getoptions(options, 'exlude_low', 1);
etgt = getoptions(options, 'etgt', []);
drawiter = getoptions(options, 'drawiter',0);
wlun = getoptions(options, 'wlun', []);

%% initial guess
if isfield(options, 'x') && not(isempty(options.x))
    x = options.x;
else
    x = applyop( U,D, y,-1,options  );
end

%% target oracle
x0 = getoptions(options, 'x0', []);
f0 = [];
if not(isempty(x0))
    f0 = applyop_single(D,x0,+1,options);
end
f0 = getoptions(options, 'f0', f0);


err = []; lun = []; Err = [];
ndisp = getoptions(options, 'ndisp', 100);
idisp = max(round(niter/ndisp),2);
verb = getoptions(options, 'verb', 1);

%% special case for BP
if T==0
    % Douglas rachford
    %	x^(t+1/2) = P_C(x^(t)) = x^(t) + pinv(A) (y - A x^(t)) 
    %	x^(t+1)   = x^(t) + mu_t (ST[2 P_C(x^(t)) - x^(t)] - x^(t+1/2)),
    gamma = getoptions(options, 'gamma', 1);    
    mu = getoptions(options, 'mu', 1);
    pA = [];
    if isnumeric(U)
        pA = pinv(U);
        if not(isempty(D))
            pA = pinv(D)*pA;
        end
    end
    
    for i=1:niter
        if verb
            progressbar(i,niter);
        end
        
        r   = y - applyop( U,D, x,+1,options, pA  );
        err(i) = norm(r, 'fro');
        
        if not(iscell(x))
            x1  = x + applyop( U,D, r,-2,options, pA  );
        else
            x1 = applyop( U,D, r,-2,options, pA  );
            for k=1:length(x)
                x1{k}  = x{k} + x1{k};
            end
        end
        
        if iscell(x1)
            lun(i) = compute_lun({x1{1:end-1}});
        else
            lun(i) = sum( abs(x1(:)) );
        end
        
        if not(iscell(x))
            x = x - mu*x1 + mu*perform_thresholding(2*x1-x,gamma,thresh);
        else
            for k=1:length(x)
                x{k} = x{k} - mu*x1{k} + mu*perform_thresholding(2*x1{k}-x{k},gamma,thresh);
            end
        end
               
        % exlucde low frequency
        if iscell(x1) && exlude_low
            x{end} = (1-mu)*x{end} + mu*x1{end};
        end
        
        % display
        f1 = applyop_single(D,x1,+1,options, pA);   
        if drawiter && mod(i,idisp)==1            
            clf; 
            if nb_dims(f1)==1
                plot( real(f1) ); 
            else
                imageplot(f1); 
            end
            drawnow;
        end
        
        if 0
        if not(isempty(f0))
            err(i) = psnr(f0,f1,1);
            % lun(i) = lunnorm({MW1{1:end-1}});
            if i==1 || err(i)>max(err(1:i-1))
                % xbest = x1;
            end            
        end
        end
        
    end
    x = x1;
    if exist('xbest')
        x = xbest;
    end
    return;
end

tau = getoptions(options, 'tau', []);
if isempty(tau)
    if isnumeric(U)
        tau = .5 * norm(U).^2;
    else not(isnumeric(U))
        warning('You should provide options.tau');
        tau = 1;
    end
end

%% MCA-like with decaying threshold
x1 = x;
for i=1:niter
    if verb
        progressbar(i,niter);
    end
    
    % residual    
    r   = y - applyop( U,D, x1,+1,options  );
    err(i) = norm(r, 'fro');
        
    % update T value
    if not(isempty(etgt)) && i>20
        mu = etgt / norm(r, 'fro' );
        mu = max(min(mu,1.1),1/1.1);
        T = T * mu;
        Tlist(i) = T;
    else
        T = Tlist(i);
    end    
    
    % x <- Thresh_{T/tau}( x + 1/tau * (U*D)'*( y - (U*D)*x ) )    
    if not(iscell(x))
        x1  = x1 + 1/tau * applyop( U,D, r,-1,options  );
    else
        xx = applyop( U,D, r,-1,options  );
        x1 = cell_add(x1,xx,1,1/tau);
    end
    xx = x1;
    if not(isempty(wlun)) && not(iscell(x1))
        x1 = perform_thresholding(x1,wlun*T/tau,thresh);
    else
        x1 = perform_thresholding(x1,T/tau,thresh);
    end
    % exclude low frequency
    if iscell(x) && exlude_low
        x1{end} = xx{end};
    end
    
    % display
    if drawiter && mod(i,idisp)==1
        f1 = applyop_single(D,x1,+1,options);
        clf;
        if nb_dims(f1)==1
            plot( real(f1) );
        else
            imageplot(f1);
        end
        drawnow;
        if not(isempty(f0))
            % f1 = applyop_single(D,x1,+1,options);
            Err(end+1) = snr(f0,f1);
        end
    end
    
    
    lun(i) = compute_lun(x1);    
end
x = x1;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = applyop(U,D,x, dir,options, pA)

if nargin<6
    pA = [];
end

if dir==1
    y = applyop_single(D,x,+1,options, pA);
    y = applyop_single(U,y,+1,options, pA);
else
    y = applyop_single(U,x,dir,options, pA);
    y = applyop_single(D,y,dir,options, pA);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = applyop_single(A,x,dir,options, pA)

if isempty(A)
    y = x;
    return;
end
if isnumeric(A)
    if dir==1
        y = A*x;
    elseif dir==-1
        y = A'*x;
    else
        if exist('pA') && not(isempty(pA))
            y = pA*x;
        else
            y = pinv(A)*x;
        end
    end
else
    y = feval( A,  x, dir, options );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function callback_identity
function lun = compute_lun(x)

if iscell(x)
    lun = 0;
    for i=1:length(x)
        lun = lun + compute_lun(x{i});
    end
    return;
end
lun = sum(abs(x(:)));

function d = nb_dims(x)

% nb_dims - debugged version of ndims.
%
%   d = nb_dims(x);
%
%   Copyright (c) 2004 Gabriel Peyré

if isempty(x)
    d = 0;
    return;
end

d = ndims(x);

if d==2 && (size(x,1)==1 || size(x,2)==1)
    d = 1;
end