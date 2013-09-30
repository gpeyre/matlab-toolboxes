function sol = SolveBP(A, y, N, maxIters, lambda, OptTol)
% SolveBP: Solves a Basis Pursuit problem
% Usage
%	sol = SolveBP(A, y, N, maxIters, lambda, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%	maxIters    maximum number of PDCO iterations to perform, default 20.
%   lambda      If 0 or omitted, Basis Pursuit is applied to the data, 
%               otherwise, Basis Pursuit Denoising is applied with 
%               parameter lambda (default 0). 
%	OptTol      Error tolerance, default 1e-3
% Outputs
%	 sol        solution of BP
% Description
%   SolveBP solves the basis pursuit problem
%      min ||x||_1 s.t. A*x = y
%   by reducing it to a linear program, and calling PDCO, a primal-dual 
%   log-barrier algorithm. Alternatively, if lambda ~= 0, it solves the
%   Basis Pursuit Denoising (BPDN) problem 
%      min lambda*||x||_1 + 1/2||y - A*x||_2^2
%   by transforming it to an SOCP, and calling PDCO.  
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as a function. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
% See Also
%   SolveLasso, SolveOMP, SolveITSP
%

if nargin < 6,
	OptTol = 1e-3;
end
if nargin < 5,
	lambda = 0;
end
if nargin < 4,
    maxIters = 20;
end

n = length(y);

n_pdco = 2*N;    % Input size
m_pdco = n;      % Output size

% upper and lower bounds
bl = zeros(n_pdco,1);
bu = Inf .* ones(n_pdco,1);

% generate the vector c
if (lambda ~= 0)
    c = lambda .* ones(n_pdco,1);
else
    c = ones(n_pdco,1);
end

% Generate an initial guess
x0 = ones(n_pdco,1)/n_pdco;       % Initial x
y0 = zeros(m_pdco,1);             % Initial y
z0 = ones(n_pdco,1)/n_pdco;       % Initial z

d1 = 1e-4;                 % Regularization parameters
if (lambda ~= 0) % BPDN
    d2 = 1; 
else
    d2 = 1e-4;
end

xsize = 1;                 % Estimate of norm(x,inf) at solution
zsize = 1;                 % Estimate of norm(z,inf) at solution

options = pdcoSet;         % Option set for the function pdco
options = pdcoSet( options, ...
                     'MaxIter    ', maxIters  , ...
                     'FeaTol     ', OptTol    , ...
                     'OptTol     ', OptTol    , ...
                     'StepTol    ', 0.99      , ...
                     'StepSame   ', 0         , ...
                     'x0min      ', 0.1       , ...
                     'z0min      ', 1.0       , ...
                     'mu0        ', 0.01      , ...
                     'method     ', 1         , ...
                     'LSQRMaxIter', 20        , ...
                     'LSQRatol1  ', 1e-3      , ...
                     'LSQRatol2  ', 1e-15     , ... 
                     'wait       ', 0    );

if (ischar(A) || isa(A, 'function_handle'))
    [xx,yy,zz,inform,PDitns,CGitns,time] = ...
        pdco(c, @pdcoMat, y, bl, bu, d1, d2, options, x0, y0, z0, xsize, zsize);
else
    Phi = [A -A];
    [xx,yy,zz,inform,PDitns,CGitns,time] = ...
        pdco(c, Phi, y, bl, bu, d1, d2, options, x0, y0, z0, xsize, zsize);
end

% Extract the solution from the output vector x
sol = xx(1:N) - xx((N+1):(2*N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    function y = pdcoMat(mode,m,n,x)
        if (mode == 1) % Direct operator
            % Decompose input
            n2 = n/2;
            u = x(1:n2);
            v = x(n2+1:n);

            % Apply matrix A
            Au = feval(A,1,m,n2,u,1:n2,n2);
            Av = feval(A,1,m,n2,v,1:n2,n2);

            y = Au-Av;
        else % Adjoint operator
            n2 = n/2;
            Atx = feval(A,2,m,n2,x,1:n2,n2);
            y = [Atx; -Atx];
        end
    end

end

%
% Copyright (c) 2006. Yaakov Tsaig
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
