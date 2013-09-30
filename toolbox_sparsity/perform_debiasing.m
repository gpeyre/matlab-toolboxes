function [x1,err] = perform_debiasing(A,x,y, options)

% perform_debiasing - remove bias by orthogonal projection
%
%   x1 = perform_debiasing(A,x,y, options);
%
%   Compute x1 with same support I=find(abs(x)>Thresh) as x that minimize
%       min | y - A(:,I)*x(I) |
%   Thresh is set in options.Thresh
%
%   Usefull to remove some bias after L1 minimization.
%
%   If A is an implicit callback function or if options.use_pinv=0, then it
%   used an itertive gradient descent.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
use_pinv = getoptions(options, 'use_pinv', 1);
Thresh = getoptions(options, 'Thresh', 1e-5);
verb = getoptions(options, 'verb', 1);

if isnumeric(A) && use_pinv
    if size(x,2)>1
        % multiple signals
        x1 = x*0;
        for i=1:size(x,2)
            [x1(:,i),err] = perform_debiasing(A,x(:,i),y(:,i), options);
        end
        return;
    end
    % use explicit pseudo inverse
    I = find(abs(x)>Thresh);
    x1 = x*0;
    x1(I) = A(:,I)\y;
    err = [];
    return;
end

niter = getoptions(options, 'niter_debiasing', 100);
tau = getoptions(options, 'tau', 1);
x1 = getoptions(options, 'xguess', x);

err = [];
x1 = perform_support_projection(x1,x,Thresh);
for i=1:niter
    if verb
        progressbar(i,niter);
    end
    % gradient descent
    % x <- x + 1/tau * A'*( y-A*x )
    r   = y - applyop_single( A, x1,+1,options  );
    err(i) = norm(r, 'fro');
    if not(iscell(x))
        x1  = x1 + 1/tau * applyop_single( A, r,-1,options  );
    else
        xx = applyop_single( A, r,-1,options  );
        for k=1:length(x)
            x1{k}  = x1{k} + 1/tau * xx{k};
        end
    end
    % project on support
    x1 = perform_support_projection(x1,x,Thresh);
end



%%
function x = perform_support_projection(x,x0,Thresh)

if iscell(x)
    for i=1:length(x)
        x{i} = perform_support_projection(x{i},x0{i},Thresh);
    end
    return;
end
% impose same support
x(abs(x0)<Thresh) = 0;

%%
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
