function rec = perform_bp_benchmark(P,W, slist, ntrial)

% perform_bp_benchmark - benchmark basis pursuit recovery
%
% rec = perform_bp_benchmark(P,W, slist, ntrial);
%
%   Try to solve many (ntrial) problems of the kind
%       x = argmin_a |a|_1    s.t. U*a=U*x0
%   where x_0 is an s-sparse signal.
%   U = P*W
%   The report the number of coordinate of x0 acuatly recovered in x.
%
%   rec(i) is the average ratio of recoverred dirac for sparsity slist(i).
%   rec=1 means perfect recovery and rec=0 means no recovery (bad).
%
%   Copyright (c) 2007 Gabriel Peyr?

if nargin<2
    slist = 10:18;
end
if nargin<3
    ntrial = 50;
end

% BP options
maxIters = 20; lambda = 0; OptTol = 1e-5;

U = P*W;
m = size(U,2);
rec = [];


% detection threshold
T = 0.1;

for is = 1:length(slist)
    progressbar(is,length(slist));
    % sparsity
    s = slist(is);
    e = 0;
    for k = 1:ntrial
        x = compute_rand_sparse(m,s);
        y = U*x;
        x0 = SolveBP(U, y, m, maxIters, lambda, OptTol);
        % compute the number of coordinates recovered
%        [tmp,I] = sort(x0);
%        x0 = x0 * 0; x0(I(end-s+1:end)) = 1;
%        e = e + sum(x0.*x)/s;
        e = e + norm( (W*x)-(W*x0) );
    end
    rec(end+1) = e/ntrial;
end