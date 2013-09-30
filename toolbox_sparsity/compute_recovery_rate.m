function rec = compute_recovery_rate(x,x0)

% compute_recovery_rate - compute the number of recovered coefficients
%
%   rec = compute_recovery_rate(x,x0);
%
%   x is the original signal with 0/1 coefficients
%   x0 is the recovered signal.
%
%   rec=1 means perfect recovery and rec=0 means no recovery (bad).
%
%   Copyright (c) 2007 Gabriel Peyre

% works for positive coefficients 0/1

s = sum(x);
x = double( x(:)>0 );

[tmp,I] = sort(x0(:));
x0 = x0 * 0; x0(I(end-s+1:end)) = 1;
rec = sum(x0.*x)/s;