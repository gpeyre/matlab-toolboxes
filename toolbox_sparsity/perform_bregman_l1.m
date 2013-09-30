function Out = perform_bregman_l1(n,A,b,mu,M,opts,varargin)


% Solve the problem
%   min ||x||_1, subject to Ax = b
% by calling the solver FPC for solving multiple instances of
%   min mu*||x||_1 + 0.5*||Ax-b^k||^2 .
% (FPC can be substituted by other solvers for the same subproblem)
%
% Technical Report: 
%   W. Yin, S. Osher, D. Goldfarb, and J. Darbon. 
%   Bregman Iterative Algorithms for l1-Minimization with Applications to Compressed Sensing.
%   Rice CAAM TR07-13 or UCLA CAM 07-37
%
% INPUT ARGUMENTS:
% n, A, b, mu: as shown above
% M: see fpc solver
%   if M=[], use 2-norm for A*x-b;
%   if M is not empty, use M-norm for A*x-b
%   (If you don't want how to use M, just let M=[])
% opts:
%   a structure of options, see l1_Bregman_opts.m
% varargin:
%   placeholder for additional parameters that are passed to FPC
%   (If you don't want how to use varargin, just ignore it)
%
% Basic Output:
%   Out.x               the last iterate x^k
%   Out.dual            an approx. dual solution, i.e., a solution for
%                       max_y <b,y>,  subject to norm(A'y,'inf') <= 1
%   Out.Bregman_itr     number of total Bregman iterations 
%                       =inf if the max # of iterations is reached before convergence
%   Out.Out_fpc         the output of last call to FPC (see fpc.m or fpc_bb.m)
%
% Copyright 2007. Wotao Yin, CAAM, Rice University
%

%% define constants
inv_norm_b = 1/norm(b);
if (any(opts.disp & [1 0 0 0 0]) || any(opts.record & [0 1 0 0 0 0])) && ~isempty(opts.fpc_opts.xs);
    inv_norm_xs = 1/norm(opts.fpc_opts.xs); 
    nz_tol=norm(opts.fpc_opts.xs,inf)*opts.rel_nz_eps; 
    nz_xs=abs(opts.fpc_opts.xs)>nz_tol; 
end
if strcmp(opts.fpc,'std'); fpc_solver = @fpc;
    elseif strcmp(opts.fpc,'bb'); fpc_solver = @fpc_bb;
    else error('no fpc solver'); end
if (opts.inv_mu); mu = 1/mu; end;
if (opts.use_rel_tol); tol = opts.rel_res_tol*norm(b); else tol = opts.abs_res_tol; end
if (opts.var_fpc_gtol); opts.fpc_opts.gtol = opts.var_fpc_gtol_max; end

% data recording initialization
if any(opts.record)
    if any(opts.record & [0 0 0 0 0 1]); Out.fpc_itrs = zeros(opts.itr_max,1); Out.fpc_time = zeros(opts.itr_max,1); end
    if any(opts.record & [0 0 0 0 1 0]); Out.fpc_res = zeros(opts.itr_max,1); Out.fpc_res_abs = zeros(opts.itr_max,1); end
    if any(opts.record & [0 0 0 1 0 0]); Out.fpc_gap = zeros(opts.itr_max,1); Out.fpc_res_err1 = zeros(opts.itr_max,1); Out.fpc_res_err2 = zeros(opts.itr_max,1); end
    if any(opts.record & [0 0 1 0 0 0]); Out.fpc_gtol = zeros(opts.itr_max,1); end
    if any(opts.record & [0 1 0 0 0 0]); Out.fpc_err = zeros(opts.itr_max,1); Out.fpc_err_abs = zeros(opts.itr_max,1); end
    if any(opts.record & [1 0 0 0 0 0]); Out.fpc_xk = []; Out.fpc_bk = []; end
end

% data dimension
m = length(b);
% implicit or explicit A
imp_A = isa(A,'function_handle');
if ~imp_A && isempty(M) && opts.check_eigmax
    eigs_opts.tol = 1E-4; eigs_opts.disp = 0; eigs_opts.issym = true;
    fprintf(1, 'checking the eigenvalue of AA''==1...');
    if abs(eigs(A*A',1,'lm',eigs_opts)-1)>1e-4; error('AA'' does not have eigenvalue 1, not satisfying the requirement by the subproblem solver FPC.'); end
    fprintf(1, 'OK!\n\n');
end

% main Bregman iterations
redo = false;
bk = b;
for k=1:opts.itr_max
    % solve the subproblem by calling fpc
    fpc_t = cputime;
    Out_fpc = fpc_solver(n,A,bk,mu,M,opts.fpc_opts,varargin{:});
    fpc_t = cputime - fpc_t;
    % data process
    if (imp_A); Axk = A(false,m,n,Out_fpc.x,[],varargin{:}); else Axk = A*Out_fpc.x; end
    b_Axk = b - Axk; bk_Axk = bk - Axk;
    res_norm = norm(b_Axk,2);
    xk_norm_1 = norm(Out_fpc.x,1);
    xk_norm_inf = norm(Out_fpc.x,inf);
    gap = abs((b'*bk_Axk)*mu - xk_norm_1);
    err1 = (b_Axk'*bk_Axk)*mu;
    err2 = (bk_Axk'*Axk)*mu - xk_norm_1;
    if (any(opts.disp & [1 0 0 0 0]) || any(opts.record & [0 1 0 0 0 0]))...
       && ~isempty(opts.fpc_opts.xs); err=norm(opts.fpc_opts.xs-Out_fpc.x); end
    % display
    if any(opts.disp)
        fprintf(1,'Itr %d:',k);
        if any(opts.disp & [0 0 0 0 1]);     fprintf(1,' itr=%3d time=%.2fs',Out_fpc.itr,fpc_t); end
        if any(opts.disp & [0 0 0 1 0]);     fprintf(1,' l1=%.2e res=%.2e (%.3e)',norm(Out_fpc.x,1),res_norm,res_norm*inv_norm_b); end
        if any(opts.disp & [0 0 1 0 0]);     fprintf(1,' gap=%.2e e1=%+.2e e2=%+.2e',gap,err1,err2); end
        if any(opts.disp & [0 1 0 0 0]);    fprintf(1,' gtol=%.2e',opts.fpc_opts.gtol); end
        if any(opts.disp & [1 0 0 0 0]) && exist('err','var')
            nz_xk = abs(Out_fpc.x)>nz_tol;
            fprintf(1,' err=%.2e (rel=%.2e) miss=%d over=%d',err,err*inv_norm_xs,nnz(nz_xs & (~nz_xk)),nnz(nz_xk & (~nz_xs))); 
        end
        fprintf(1,'.\n\n');
    end
    % data recording
    if any(opts.record)
        if any(opts.record & [0 0 0 0 0 1]); Out.fpc_itrs(k) = Out_fpc.itr; Out.fpc_time(k) = fpc_t; end
        if any(opts.record & [0 0 0 0 1 0]); Out.fpc_res(k) = res_norm; Out.fpc_res_abs(k) = res_norm*inv_norm_b; end
        if any(opts.record & [0 0 0 1 0 0]); Out.fpc_gap(k) = gap; Out.fpc_res_err1(k) = err1; Out.fpc_res_err2(k) = err2; end
        if any(opts.record & [0 0 1 0 0 0]); Out.fpc_gtol(k) = opts.fpc_opts.gtol; end
        if any(opts.record & [0 1 0 0 0 0]) && exist('err','var'); Out.fpc_err(k) = err; Out.fpc_err_abs(k) = err*inv_norm_xs; end
        if any(opts.record & [1 0 0 0 0 0]); Out.fpc_xk = [Out.fpc_xk Out_fpc.x]; Out.fpc_bk = [Out.fpc_bk bk]; end
    end
    % convergence test
    if (res_norm<tol); break; end
    % failure test
    if (Out_fpc.itr == inf); disp('FPC fail to converge; please use larger mu or allow more FPC iterations.'); break; end
    if (Out_fpc.itr == 0); disp('mu is too large; FPC returns the solution 0; please reduce mu.'); break; end
    % data update
    if abs(err2) <= abs(err1)*5;
        % proceed to next iteration
        redo = false;
        bk = b + bk_Axk;
        if (opts.use_last); opts.fpc_opts.x0 = Out_fpc.x; end
        if (opts.var_fpc_gtol);
            opts.fpc_opts.gtol = max(min( abs(err1/mu)/(nnz(Out_fpc.x)*xk_norm_inf), opts.var_fpc_gtol_max),opts.var_fpc_gtol_min); 
        end
    else
        % need to redo last step
        if redo; disp('Fail to decrease error. Exiting...'); break; end
        redo = true; disp('Rerun last step with smaller gtol.');
        opts.fpc_opts.gtol = opts.fpc_opts.gtol*abs(err1/err2);
        if Out_fpc.itr == inf; opts.fpc_opts.mxitr = min(round(opts.fpc_opts.mxitr*abs(err2/err1)),opts.fpc_max_mxitr); end
    end
end

% save output
if (res_norm<tol); Out.Bregman_itr = k; else Out.Bregman_itr = inf; end
Out.dual = bk_Axk*mu;
Out.Out_fpc = Out_fpc;
Out.x = Out_fpc.x;

end



% Options for l1 Bregman iterations
% opts:                     
%   opts.fpc                              % fpc solver selection: 'std', 'bb'
%   opts.fpc_opts                         % options passed to fpc
%   opts.check_eigmax                     % checking max eigenvale of AA' == 1 or nor (to satisfy FPC's requirement)
%   opts.fpc_max_mxitr                    % the ceiling of mxitr for fpc
%   opts.itr_max                          % maximum Bregman iterations. Default: 50
%   opts.use_rel_tol                      % true: use opts.rel_res_tol; false; use abs_res_tol. Default: true
%     opts.rel_res_tol                    % tolerance of relative residual in 2-norm, ||Au-f||/||f||
%     opts.abs_res_tol                    % tolerance of absolute residual in 2-norm, ||Au-f||
%   opts.inv_mu                           % need to invert mu for the subproblem solver
%   opts.var_fpc_gtol                     % varying gtol. Default: true
%     opts.var_fpc_gtol_max               % max gtol if var_fpc_gtol == true. Default: 1e-2
%     opts.var_fpc_gtol_min               % min gtol if var_fpc_gtol == true. Default: 1e-10
%   opts.use_last                         % use last solution as next starting point, suggest: false (not improvement observed with "true")
%   opts.rel_nz_eps                       % used to determine zeros in x^k relative to norm(x^k,inf). Default: 1e-3;
%   opts.disp                             % display options (see below). Default: no display
%   opts.record                           % data recording options (see below). Default: bin 000001
%
%
% Display options:
%   bit 0 right - basic information: fpc_itr, fpc_time
%   bit 1 - norm(x,1) and residual info: norm(A*x^k-b), norm(A*x^k-b)/norm(b)
%   bit 2 - accuracy info: gap abs(||x^k||_1 - <b^k-A*x^k,b>/mu), err1 <b^k-A*x^k,b-A*x^k>/mu, err2 <b^k-A*x^k,A*x^k>/mu - ||x^k||_1
%   bit 3 - gtol
%   bit 4 - error. if opts.fpc.xs exists, norm(xs-x^k), norm(xs-x^k)/norm(xs), missed nz, over nz
% Data recording options:
%   bit 0 right - basic information: fpc_itr, fpc_time
%   bit 1 - residual info
%   bit 2 - accuracy info
%   bit 3 - gtol
%   bit 4 - error
%   bit 5 - A, b, x^k, b^k
%
%  Copyright 2007. Wotao Yin, CAAM, Rice University

function opts = l1_Bregman_opts(opts)

%% default values
fpc_default = 'bb';
fpc_max_mxitr_default = 1e4;
itr_max_default = 50;
check_eigmax_default = true;
use_rel_tol_default = true;
    rel_res_tol_default = 1e-5;
    abs_res_tol_default = 1e-5;
inv_mu_default = true;
var_fpc_gtol_default = true;
    var_fpc_gtol_max_default = 1e-2;
    var_fpc_gtol_min_default = 1e-10;
    fix_fpc_gtol_default = 1e-5;
use_last_default = false;
rel_nz_eps_default = 1e-3;
disp_default = [0 0 0 0 0];
record_default = [0 0 0 0 0 1];

err_msg=[];

%% option assignment
if ~exist('opts','var'); opts = []; end

if ~isfield(opts,'fpc'); opts.fpc = fpc_default; end
if ~isfield(opts,'fpc_opts') && strcmp(opts.fpc,'std'); opts.fpc_opts = fpc_opts([]); end
if ~isfield(opts,'fpc_opts') && strcmp(opts.fpc,'bb'); opts.fpc_opts = fpc_bb_opts([]); end

if ~isfield(opts,'check_eigmax'); opts.check_eigmax = check_eigmax_default; end
if ~isfield(opts,'fpc_max_mxitr'); opts.fpc_max_mxitr = fpc_max_mxitr_default; end
if ~isfield(opts,'itr_max'); opts.itr_max = itr_max_default; end 
if ~isfield(opts,'opts.use_rel_tol'); opts.use_rel_tol = use_rel_tol_default; end
if opts.use_rel_tol == true && ~isfield(opts,'rel_res_tol'); opts.rel_res_tol=rel_res_tol_default; end
if opts.use_rel_tol == false && ~isfield(opts,'abs_res_tol'); opts.abs_res_tol=abs_res_tol_default; end
if ~isfield(opts,'inv_mu'); opts.inv_mu = inv_mu_default; end
if ~isfield(opts,'var_fpc_gtol'); opts.var_fpc_gtol = var_fpc_gtol_default; end
    if (opts.var_fpc_gtol == true && ~isfield(opts,'var_fpc_gtol_max')); opts.var_fpc_gtol_max=var_fpc_gtol_max_default; end
    if (opts.var_fpc_gtol == true && ~isfield(opts,'var_fpc_gtol_min')); opts.var_fpc_gtol_min=var_fpc_gtol_min_default; end
    if opts.var_fpc_gtol == false
        if ~isfield(opts,'fix_fpc_gtol'); opts.fpc_opts.gtol = fix_fpc_gtol_default;
        else opts.fpc_opts.gtol = opts.fix_fpc_gtol; end
    end
if ~isfield(opts,'use_last'); opts.use_last = use_last_default; end
if ~isfield(opts,'rel_nz_eps'); opts.rel_nz_eps = rel_nz_eps_default; end
if ~isfield(opts,'disp'); opts.disp = disp_default; end
if ~isfield(opts,'record'); opts.record = record_default; end

% if ~isempty(err_msg)
%     disp([err_msg(3:end) '.']);
%     error('Error in opts.');
% end


% Options for Fixed Point Continuation (FPC)
%
%--------------------------------------------------------------------------
% DESCRIPTION
%--------------------------------------------------------------------------
%
% opts = fpc_opts(opts)
%
% If opts is empty upon input, opts will be returned containing the default
% options for fpc.m.  
%
% Alternatively, if opts is passed with some fields already defined, those
% fields will be checked for errors, and the remaining fields will be added
% and initialized to their default values.
%
% Table of Options.  ** indicates default value.
%
% FIELD   OPTIONAL  DESCRIPTION
% .x0       YES     Initial value of x.  If not defined, x will be
%                   initialized according to .init.
% .xs       YES     Original signal xs.  If passed, fpc will calculate and
%                   output ||x - xs||/||xs|| in vector Out.n2re.
% .init     YES     If .x0 is not defined, .init specifies how x is to be
%                   initialized.  
%                       0 -> zeros(n,1)
%                       1 -> x = tau*||AtMb||_Inf * ones(n,1)
%                    ** 2 -> x = tau*AtMb **
% .tau      YES     0 < .tau < 2.  If not specified, tau is initialized 
%                   using a piecewise linear function of delta = m/n.
% .mxitr    NO      Maximum number of inner iterations.
%                   ** 1000 **
% .eta      NO      Ratio of current ||b - Ax|| to approx. optimal 
%                   ||b - Ax|| for next mu value.
%                   ** 4 **
% .fullMu   NO      If true, then mu = eta*sqrt(n*kap)/phi, where phi = 
%                   ||Ax - b||_M, which guarantees that phi(now)/phi(next) 
%                   >= eta.  Otherwise mu = eta*mu.
%                   ** false **
% .kappa    YES     Required if fullMu.  Is ratio of max and min 
%                   eigenvalues of M^{1/2}AA'M^{1/2} (before scaling).
%                   ** not supplied **
% .xtol     NO      Tolerance on norm(x - xp)/norm(xp).
%                   ** 1E-4 **
% .gtol     NO      Tolerance on mu*norm(g,'inf') - 1
%                   ** 0.2 **
%--------------------------------------------------------------------------

function opts = fpc_opts(opts)

if isfield(opts,'x0')
    if ~isvector(opts.x0) || ~min(isfinite(opts.x0))
        error('If used, opts.x0 should be an n x 1 vector of finite numbers.');
    end
elseif isfield(opts,'init')
    if (opts.init < 0) || (opts.init > 2) || opts.init ~= floor(opts.init)
        error('opts.init must be an integer between 0 and 2, inclusive.');
    end
else
    opts.init = 2; 
end

if isfield(opts,'xs')
     if ~isvector(opts.xs) || ~min(isfinite(opts.xs))
        error('If passed, opts.xs should be an n x 1 vector of finite numbers.');
     end
end

if isfield(opts,'tau')
    if (opts.tau <= 0) || (opts.tau >= 2)
        error('If used, opts.tau must be in (0,2).');
    end
end
    
if isfield(opts,'mxitr')
    if opts.mxitr < 1 || opts.mxitr ~= floor(opts.mxitr)
        error('opts.mxitr should be a positive integer.');
    end
else
    opts.mxitr = 1000;
end
    
if isfield(opts,'xtol')
    if (opts.xtol < 0) || (opts.xtol > 1)
        error('opts.xtol is tolerance on norm(x - xp)/norm(xp).  Should be in (0,1).');
    end
else
    opts.xtol = 1E-4;
end

if isfield(opts,'gtol')
    if (opts.gtol < 0) || (opts.gtol > 1)
        error('opts.gtol is tolerance on mu*norm(g,''inf'') - 1.  Should be in (0,1).');
    end
else
    opts.gtol = 0.2;
end

if isfield(opts,'eta')
    if opts.eta <= 1
        error('opts.eta must be greater than one.');
    end
else
    opts.eta = 4;
end

if isfield(opts,'fullMu')
    if ~islogical(opts.fullMu)
        error('fullMu should be true or false.');
    end
else
    opts.fullMu = false;
end

if isfield(opts,'kappa')
    if opts.kappa < 1
        error('opts.kappa is a condition number and so should be >= 1.');
    end
end

return

% Copyright (c) 2007.  Elaine Hale, Wotao Yin, and Yin Zhang
%
% Last modified 28 August 2007.


% Fixed Point Continuation (FPC) for l1 Regularized Least Squares
%
%--------------------------------------------------------------------------
% GENERAL DESCRIPTION & INPUTS
%-------------------------------------------------------------------------- 
%
% Out = fpc(n,A,b,mu,M,opts,varargin)
%
% Solves
%
%   min ||x||_1 + (mu/2)*||Ax - b||_M^2.
%
% A may be an explicit m x n matrix, or a function handle that implements 
% A*x and A'*x.  If the latter, this function must have the form
%
%   function y = name(trans,m,n,x,inds,varargin)
%
% where
%
%   trans    - if false then y = A(:,inds)*x, if true then y = A(:,inds)'*x 
%   m, n     - A is m x n
%   x        - length(inds) x 1 vector if ~trans; m x 1 vector if trans
%   inds     - vector of indices; if empty the function should return A*x 
%              or A'*x, as appropriate
%   varargin - placeholder for additional parameters
%
% b must be an m x 1 vector.
%
% M can be any m x m positive definite matrix, or the empty matrix.  If M
% is empty, fpc assumes M = I, which reduces the second term of the 
% objective to (mu/2)*||Ax - b||_2^2 (standard least squares). 
%
% This function assumes that the maximum eigenvalue of A^T M A is less
% than or equal to 1.  If your initial problem does not satisfy this
% condition, an equivalent problem may be constructed by setting:
%   sigma^2 = max eigenvalue of A^T M A
%   mu = mu*sigma^2
%   A = A/sigma
%   b = b/sigma
% Also see getM_mu.m.
%
% fpc_opts.m describes the available options.  If opts is passed as empty,
% fpc_opts will be called to obtain the default values.
%
% All variables in varargin are passed to A if A is a function handle.
%
%--------------------------------------------------------------------------
% OUTPUTS
%--------------------------------------------------------------------------
%
%   Out.x    - x at last iteration
%   Out.f    - vector of function values
%   Out.lam  - vector of ||x||_1
%   Out.step - vector of norm(x - xp)
%   Out.mus  - vector of mu values for each outer iteration
%   Out.itr  - number of iterations to convergence (or Inf if reach mxitr)
%   Out.itrs - vector of number of inner iterations completed during each
%              outer iteration
%   Out.tau  - value of tau
%   Out.n2re - if opts.xs exists, is vector of norm(x - xs)/norm(xs).
%              starts with 0th iteration.
%--------------------------------------------------------------------------

function Out = fpc(n,A,b,mu,M,opts,varargin)

% problem dimension
m = length(b);

% implicit or explicit A
imp = isa(A,'function_handle');
% calculate AtMb
if imp
    if isempty(M), AtMb = A(true,m,n,b,[],varargin{:});
    else AtMb = A(true,m,n,M*b,[],varargin{:}); end
else
    if isempty(M), AtMb = A'*b;
    else AtMb = A'*(M*b); end
end

% check for 0 solution
if mu <= 1/norm(AtMb,'inf'); 
    Out.x = zeros(n,1); Out.itr = 0; Out.itrs = 0; 
    Out.tau = 0; Out.mus = mu; Out.lam = 0; Out.step = [];
    if isempty(M), Out.f = (mu/2)*(b'*b); 
    else Out.f = (mu/2)*(b'*(M*b)); end
    if isfield(opts,'xs'), Out.n2re = 1; end
    return
end

% get opts
if isempty(opts), opts = fpc_opts([]); end

% initialize x, nu, tau, mu
muf = mu;                       % final value of mu
[x,nu,tau,mu] = fpc_init(n,m,b,AtMb,M,opts);
if mu > muf, mu = muf; nu = tau/mu; end
Out.mus = mu; Out.tau = tau;

% initialize Out.n2re
if isfield(opts,'xs'), xs = opts.xs; else xs = []; end
if ~isempty(xs), Out.n2re = norm(x - xs)/norm(xs); end

xtol = opts.xtol;
gtol = opts.gtol;

% prepare for iterations
Out.step = []; Out.itrs = []; Out.f = []; Out.lam = [];
Out.itr = Inf; oitr = 0; Ax = [];

% main loop
for i = 1:opts.mxitr
    
    % store old point
    xp = x;
    
    % get gradient at x and store objective function
	[g,f,lam] = get_g(x,imp,m,n,mu,A,b,M,AtMb,varargin{:});
    Out.f = [Out.f; f]; Out.lam = [Out.lam; lam];
    
    % take fixed-point step
    y = x - tau*g; 
    x = sign(y).*max(0,abs(y)-nu);
    
    nrmxxp = norm(x - xp);
    Out.step = [Out.step; nrmxxp]; 
    
    if ~isempty(xs), Out.n2re = [Out.n2re; norm(x - xs)/norm(xs)]; end
    
    crit1 = nrmxxp/max(norm(xp),1);
    crit2 = mu*norm(g,'inf') - 1;
    
    if (crit1 < xtol*sqrt(muf/mu)) && (crit2 < gtol)
        
        oitr = oitr + 1;
        
        if isempty(Out.itrs), Out.itrs = i;
        else Out.itrs = [Out.itrs; i - sum(Out.itrs)]; end
        
        % stop if reached muf
        if mu == muf
            Out.x = x; Out.itr = i;
            [g,f,lam] = get_g(x,imp,m,n,mu,A,b,M,AtMb,varargin{:});
            Out.f = [Out.f; f]; Out.lam = [Out.lam; lam];
            return 
        end
        
        % update mu
        if opts.fullMu
            phi = sqrt((2/mu)*(f - lam));
            mu = getNextMu(n,phi,opts.eta,opts.kappa);
        else
            mu = opts.eta*mu;
        end
        mu = min(mu,muf); nu = tau/mu;          
        Out.mus = [Out.mus; mu];
    end
end

% did not converge within opts.mxitr
Out.x = x;
if isempty(Out.itrs), Out.itrs = i;
else Out.itrs = [Out.itrs; i - sum(Out.itrs)]; end

end % fpc

%--------------------------------------------------------------------------
% SUBFUNCTION FOR INITIALIZATION
%--------------------------------------------------------------------------
%
% OUTPUTS -----------------------------------------------------------------
% x   - initialized based on opts.x0 and opts.init.  if opts.x0 exists, 
%       x = opts.x0, otherwise, opts.init determines x:
%           0 - x = zeros(n,1)
%           1 - x = tau*||AtMb||_Inf * ones(n,1)
%           2 - x = tau*AtMb
% nu  - equals tau/mu
% tau - equals opts.tau if exists, otherwise min(1.999,-1.665*m/n+2.665)
% mu  - set based on x = 0, mu = 1/norm(AtMb,inf), and getNextMu
%--------------------------------------------------------------------------

function [x,nu,tau,mu] = fpc_init(n,m,b,AtMb,M,opts)

% initialize tau
if isfield(opts,'tau'), tau = opts.tau;
else tau = min(1.999,-1.665*m/n + 2.665); end

% initialize x
if isfield(opts,'x0')
    x = opts.x0;
    if length(x) ~= n, error('User supplied x0 is wrong size.'); end
else
    switch opts.init
        case 0, x = zeros(n,1);
        case 1, x = tau*norm(AtMb,inf)*ones(n,1);
        case 2, x = tau*AtMb;
    end
end

% initialize mu
if opts.fullMu && isfield(opts,'kappa')
    if isempty(M), phi = norm(b);       % phi = ||Ax - b||_M
    else phi = sqrt(b'*(M*b)); end
    mu = getNextMu(n,phi,opts.eta,opts.kappa);
else
    if opts.fullMu
        warning('opts.fullMu = true, but opts.kappa is not supplied.  Switching to opts.fullMu = false.');
        opts.fullMu = false;
    end
    mu = opts.eta/norm(AtMb,inf);
end

% initialize nu
nu = tau/mu;

end % fpc_init

%--------------------------------------------------------------------------
% SUBFUNCTION FOR CALCULATING NEXT mu
%--------------------------------------------------------------------------
%
% Calculates the next value of mu based on taking a predictor step along
% the pareto curve phi(lam).  The derivative of this curve is derived in 
% 
% van den Berg, E. and M. Friedlander.  In pursuit of a root.  Preprint,
%   2007.
%
% The steplength is chosen so that phi(current)/phi(next) \approx eta.
% Mu is chosen so that the true phi(next) is guaranteed to be <= the 
% predicted value.
%
% INPUTS ------------------------------------------------------------------
% n     - length of x
% phi   - ||Ax - b||_M
% g     - A'M(Ax - b)
% eta   - parameter for choosing how much phi should decrease during next
%         outer iteration.  getNextMu is approximately the same as choosing
%         mu = eta*mu.
% kap   - condition number of A'MA.  See getM_mu.
%--------------------------------------------------------------------------

function mu = getNextMu(n,phi,eta,kap)

mu = eta*sqrt(n*kap)/phi;

end % getNextMu

%--------------------------------------------------------------------------
% SUBFUNCTION FOR CALCULATING g
%--------------------------------------------------------------------------

function [g,f,lam] = get_g(x,imp,m,n,mu,A,b,M,AtMb,varargin)

% get A*x
if imp
	Ax = A(false,m,n,x,[],varargin{:});
else
    Ax = A*x;
end

% calc g
if imp
    if isempty(M)
        g = A(true,m,n,Ax,[],varargin{:})-AtMb;
    else
        g = A(true,m,n,M*Ax,[],varargin{:})-AtMb;
    end
elseif isempty(M)
    g = A'*Ax - AtMb;
else
    g = A'*(M*Ax) - AtMb;
end

% calc f
r = Ax - b; lam = sum(abs(x));
if isempty(M)
    f = 0.5*mu*norm(r)^2 + lam;
else
    f = 0.5*mu*r'*(M*r) + lam;
end

end % get_g

% Copyright (c) 2007.  Elaine Hale, Wotao Yin, and Yin Zhang
%
% Last modified 28 August 2007.


% RECOMMENDED M AND mu FOR l1-REGULARIZED WEIGHTED LEAST SQUARES
%
%--------------------------------------------------------------------------
% DESCRIPTION AND INPUTS
%--------------------------------------------------------------------------
% 
% [M,mu,A,b,sig,kap,tau,M12] = getM_mu(full,mu,m,n,Ameth,A,b,sig1,sig2,alpha)
%
% Constructs M and mu as described in Hale, Yin and Zhang 2007, based
% on the noise estimates sig1 and sig2, and the statistical parameter
% alpha.  
%
%  full  - only relevant if sig1 > 0 and AA' is not a multiple of I.  In 
%          this case, and if full == true, then the full estimate for M, 
%          M = (sig1^2 AA' + sig2^2 I)^{-1} is calculated.  Otherwise M is 
%          estimated to be (sig1^2 lam_max(AA') + sig2^2)^{-1} I.  The 
%          constant is then moved to mu so that M is returned as [].
%
%  mu    - if empty, then recommended mu value is calculated.
%
%  m,n   - A is an m x n matrix.
%
%  Ameth - indicates what type of A matrix is being passed.  values >= 0
%          correspond to indices from getData.m.
%               -1 - generic A matrix
%                0 - randn(m,n) (uses analytic bound on min and max
%                    eigenvalues of AA')
%                1 - 0 and then columns are scaled to unit norm
%                2 - 0 and then QR factorization to obtain orthonormal rows
%                    (AA' = I)
%                3 - bernoulli +/- 1 distribution
%                4 - partial Hadamard matrix (AA' = nI)
%                5 - partial fourier matrix (implicit,AA'=nI,ifft = A'/n)
%                6 - partial discrete cosine matrix (implicit,AA'=I)
%                7 - partial 2-d fourier matrix (as 5,but AA'=n^2*I)
%                8 - partial 2-d discrete cosine matrix (as 6)
%          Codes -1, 1, and 3 are equivalent; 2, 6 and 8 are equivalent; 4
%          and 5 are equivalent.
%
%  A,b   - problem data.  Only accessed if Ameth <= 4 (A is explicit).
%
%  sig1, sig2 - noise level estimates.  In particular we assume that 
%          b = A*(xs + epsilon1) + epsilon2, where epsiloni is iid Gaussian
%          noise with standard deviation sigi.
%
%  alpha - mu estimate uses 1-alpha chi^2 critial value.  mu is not very
%          sensitive to alpha and alternative is provided in case the 
%          statistics toolbox is not available.  Results in Hale, Yin and 
%          Zhang 2007 used alpha = 0.5.
%
%--------------------------------------------------------------------------
% OUTPUTS
%--------------------------------------------------------------------------
%
% M     - recommended value of M.  Returned as [] if M = I.
%
% mu    - if input mu is the empty matrix, output mu is the recommended 
%         value.  otherwise, mu is untouched.
%
% A,b   - problem data.  Untouched if AA' = I or A is implicit.  Otherwise, 
%         A and b are scaled so that the maximum eigenvalue of A'*M*A = 1.
%
% sig   - approximate noise level (standard deviation) of least squares 
%         solution.  equal to sqrt(sig1^2 + sig2^2/(min eig of AA')).
%
% kap   - condition number of M12*AA'*M12.
%
% tau   - if AtMA is not well conditioned, tau is returned as 2-eps.  
%         otherwise, tau is returned empty to signal that fpc default 
%         should be used.
%
% M12   - equal to M^(1/2) if ~isempty(M).
%--------------------------------------------------------------------------

function [M,mu,A,b,sig,kap,tau,M12] = getM_mu(full,mu,m,n,Ameth,A,b,sig1,sig2,alpha)

tau = []; M12 = [];

% options for eigs
opts.tol = 1E-4;
opts.disp = 0;
opts.issym = true;

% calculate M, mu, sig, kap, tau, M12
if full && sig1 > 0 && ismember(Ameth,[-1 0 1 3])
    % most general case--must calculate full M matrix
    AAt = A*A';
    M = sig1^2*AAt; M(1:m+1:m^2) = M(1:m+1:m^2) + sig2^2;
    % invert to get M
    lsopts.SYM = true; lsopts.POSDEF = true;
    M = linsolve(M,eye(m),lsopts);
    % need matrix square root
    M12 = sqrtm(M);
    M12AAtM12 = M12*AAt*M12;
    % get eigenvalues
    smax = sqrt(eigs(M12AAtM12,1,'lm',opts));
    smin = sqrt(eigs(M12AAtM12,1,0,opts));
    if Ameth == 0
        sig = sqrt(sig1^2 + sig2^2/max((1 - sqrt(m/n))^2*n,eps^2));
    else
        sig = sqrt(sig1^2 + sig2^2/max(eigs(AAt,1,0,opts),eps^2));
    end
    kap = smax^2/smin^2;
    % no tau because full M means well-conditioned AtMA
    if isempty(mu)
        if strcmp(which('chi2inv'),'');
            % avoids chi2inv, is good approximation for large n
            mu = smax*sqrt(n*kap/m);
        else
            mu = smax*sqrt(n*kap/chi2inv(1-alpha,m));
        end
    end
else
    % M is a multiple of I
    M = [];
    switch Ameth
        case {-1,1,3}
            % unknown eigenvalues
            AAt = A*A';
            smax = sqrt(eigs(AAt,1,'lm',opts));
            smin = sqrt(eigs(AAt,1,0,opts));
            muSig = sqrt(sig1^2*smax^2 + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/smin^2);
            kap = smax^2/smin^2;
            tau = 2-eps;
        case 0
            % have tight bounds
            delta = m/n;
            smax = (1 + sqrt(delta))*sqrt(n);
            smin = max((1 - sqrt(delta))*sqrt(n),eps);
            muSig = sqrt(sig1^2*smax^2 + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/smin^2);
            kap = smax^2/smin^2;
            tau = 2-eps;
        case {2,6,8}
            % AAt = I
            smin = 1; smax = 1;
            sig = sqrt(sig1^2 + sig2^2);
            muSig = sig;
            kap = 1;
        case {4,5}
            % AAt = nI
            smin = sqrt(n); smax = sqrt(n); 
            muSig = sqrt(sig1^2*n + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/n);
            kap = 1;
            % if Ameth = 5, b was already scaled by getData since A scaling
            % is implemented in pfft_wrap_part
        case 7
            smin = n; smax = n; 
            muSig = sqrt(sig1^2*n^2 + sig2^2);
            sig = sqrt(sig1^2 + sig2^2/n^2);
            kap = 1;
            % b was already scaled by getData since A scaling is 
            % implemented in pfft2_wrap_part
    end
    % calculate recommended mu if required
    if isempty(mu)
        if strcmp(which('chi2inv'),'');
            % avoids chi2inv, is good approximation for large n
            mu = (smax/muSig)*sqrt(n*kap/m);
        else
            mu = (smax/muSig)*sqrt(n*kap/chi2inv(1-alpha,m));
        end
    end
end

% scale A and b if required
if ismember(Ameth,[-1,0,1,3,4])
    A = (1/smax)*A; b = (1/smax)*b;
end

return

% Copyright (c) 2007.  Elaine Hale, Wotao Yin, and Yin Zhang
%
% Last modified 28 August 2007.