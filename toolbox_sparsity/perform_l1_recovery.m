function x = perform_l1_recovery(D,y,options)

% perform_l1_recovery - run an L1 solver
%
%   x = perform_l1_recovery(D,y,options);
%
%   options.method can be 
%       'bp', 'omp', 'itthresh' (nb.iter is options.niter_inversion).
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

if isfield(options, 'method')
    method = options.method;
else
    method = 'bp';
end
if isfield(options, 's')
    s = options.s;
else
    s = 20;
end
if isfield(options, 'niter_inversion')
    mcaIters = options.niter_inversion;
else
    mcaIters = 200;
end


m = size(D,2);

% OMP options
maxItersOMP = s;
solFreq = 0; verbose = 0; lambdaStop = 0;

% BP options
maxIters = getoptions(options, 'niterbp', 20); 
lambda_bp = getoptions(options, 'lambda_bp', 0); 
OptTol = getoptions(options, 'tol', 1e-5);

switch lower(method)
    case 'bp'
        x = SolveBP(D, y, m, maxIters, lambda_bp, OptTol);
    case 'omp'
        % x = SolveOMP(D, y, m, maxItersOMP, lambdaStop, solFreq, verbose, OptTol);
        options.nbr_max_atoms = s;
        x = perform_omp(D,y,options);
    case 'itthresh'
        options.niter = mcaIters;
        [x,E] = perform_iterative_thresholding(D,y,options);
    otherwise
        error('Unknown method');
end