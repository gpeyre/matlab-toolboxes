function options = pdcoSet(varargin)
%pdcoSet creates or alters an options structure for pdco.m.
%It is modeled after MATLAB's original version of optimset.
%
%   options = pdcoSet (with no input arguments) creates a structure
%   with all fields set to their default values.  Each field is an
%   option (also called a parameter).
%
%   pdcoSet (with no input or output arguments) displays all options
%   and their default values.
%
%   options = pdcoSet('PARAM1',VALUE1,'PARAM2',VALUE2,...);
%   creates a structure with the specified named parameters and values.
%   Any unspecified parameters are set to [].
%   It is sufficient to use only the leading characters that uniquely
%   identify the parameter.  Case is ignored for parameter names.
%
%   NOTE: For values that are strings, correct case and the complete string
%   are required; if an invalid string is provided, the default is used.
%
%   options = pdcoSet(oldopts,'PARAM1',VALUE1,...);
%   creates a copy of oldopts with the named parameters reset to the
%   specified values.
%
%   options = pdcoSet(oldopts,newopts);
%   combines an existing structure oldopts with a new structure newopts.
%   Any parameters in newopts with non-empty values overwrite the
%   corresponding old parameters in oldopts.
%
% options.MaxIter     Maximum iterations of the primal-dual barrier method.
%                     Most problems should solve within 30 PDitns.
% options.FeaTol      Accuracy for satisfying Ax + D2 r = b, A'y + z = gobj
%                     and x - x1 = bl, x + x2 = bu, where x1, x2 > 0.
%                     1e-6 is typically small enough.
%                     1e-5 may be acceptable also.
% options.OptTol      Accuracy for satisfying x1.*z1 = 0, x2.*z2 = 0,
%                     where z = z1 - z2 and z1, z2 > 0.
%                     Typically the same as Featol.
% options.StepTol     (between 0 and 1): Controls how close each an
%                     x or z may be to reaching a bound at each step.
%                     For safety, should not be bigger than 0.99 (say)
%                     for nonlinear problems.
% options.StepSame    1 (true) if stepx and stepz should be the same
%                     (gives true Newton method for nonlinear problems);
%                     0 (false) if stepx and stepz may be different
%                     (gives better performance for linear programs).
% options.x0min       Min distance between x0 and bl or bu  AFTER SCALING.
%                     1.0 is about right for cold starts.
%                     0.1 or 0.01 may be ok if x0, z0 are "good".
% options.z0min       Min distance between abs(z0) and zero AFTER SCALING,
%                     when z0 is used to initialize z1 > 0 and z2 > 0.
%                     Typically the same as x0min.
% options.mu0         Initial mu (ABSOLUTE VALUE) for solving scaled problem.
% options.Method      Specifies how each search direction (dx,dy,dz1,dz2)
%                     should be computed.  Several methods exist for
%                     experimental purposes.  If A has fewer rows than columns
%                     (m < n) we usually solve for dy first (most of the work)
%                     and then get the other components cheaply.
%
%  Method  Solve for  Using                                         Implemented?
%     1       dy      Sparse Cholesky on (A D^2 A' + D2^2 I).           Yes
%     2       dy      Sparse QR on corresponding least-squares problem  Yes
%     3       dy      LSQR on least-squares problem                     Yes
%
%    11       dx      Sparse Cholesky on (D A'A D  + D2^2 I).           No
%    12       dx      Sparse QR on corresponding least-squares problem  No
%    13       dx      LSQR on least-squares problem                     No
%
%    21    dx,dy      Sparse LU on 2x2 KKT-type system                  No
%    23    dx,dy      SYMMLQ on same system                             No
%
%    31    dx,dy,dz1  Sparse LU on 3x3 system (only for problems with   No
%                     with vanilla bounds: 0 < x < inf).
%                     This is a HUGE system, but it is relatively
%                     well-conditioned compared to the other systems.
%
%    41 dx,dy,dz1,dz2 Sparse LU on 4x4 system with general bounds.      Yes
%                     This is an even HUGER system.
%
%    If A is an explicit sparse matrix, all methods are applicable.
%     1 is usually best (e.g. for LPs).
%     2 may be more reliable; it's there for checking.
%     3 is sometimes more efficient (e.g. for entropy problems).
%       Diagonal preconditioning is possible.
%
%    If A is an operator,
%     3 must be used.  Diagonal preconditioning is not possible.
%
%    Notes:
%    Method =  1      On ENTROPY.big, symamd can find an ordering
%                     but chol never finishes.
%    Method = 41      On ENTROPY.big, symamd never finishes.
%
% The following options control LSQR when Method = 3:
%
% options.LSQRMaxIter * min(m,n) is the maximum LSQR (CG) iterations.
% options.LSQRatol1   is the starting value of the LSQR accuracy 
%                     tolerance "atol" (if LSmethod = 3).
%                     1e-3 or 1e-4 sometimes works.
%                     1e-8 may be needed for LPs.
%                     In general, if max(Pinf,Dinf,Cinf) doesn't decrease
%                     every iteration, set atol1 to a smaller value.
% options.LSQRatol2   is the smallest value atol is reduced to.
% options.LSQRconlim  shuts LSQR down early if its matrix is ill-conditioned.
%
% options.wait = 0    means pdco should proceed to solve the problem.
%              = 1    means pdco should pause to allow interactive resetting
%                     of some of the parameters.


% pdcoSet.m is derived from optimset.m (Revision 1.14, 1998/08/17)
% in the Optimization Toolbox of The MathWorks, Inc.
%
% 28 Sep 2000: First version of pdscoSet.m.
% 30 Sep 2002: First version of pdcoSet.m derived from pdscoSet.m.
% 04 Oct 2002: StepSame introduced.
% 10 Aug 2003: mu0 is now an absolute value -- the initial mu.
% 19 Sep 2003: Method = 41 implemented.
% 22 Sep 2003: LSproblem and LSmethod replaced by Method.
%              
%
% Michael Saunders, SOL, Stanford University.


if (nargin == 0)        % Set default options.
    defoptions.MaxIter      =    30;
    defoptions.FeaTol       =  1e-6;
    defoptions.OptTol       =  1e-6;
    defoptions.StepTol      =  0.99;
    defoptions.StepSame     =     1;  % 1 for stepx == stepz (NLPs)
    defoptions.x0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
    defoptions.z0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
    defoptions.mu0          =  1e-1;  % < 1.0 better than 1.0?
    defoptions.Method       =     3;  % 3 = computed dy using LSQR
    defoptions.LSQRMaxIter  =  10.0;
    defoptions.LSQRatol1    = 1e-08;
    defoptions.LSQRatol2    = 1e-15;  % 
    defoptions.LSQRconlim   = 1e+12;  % Somewhere between e+8 and e+16
    defoptions.wait         =     0;
    defoptions.NOTE         = 'LSQRMaxIter is scaled by the matrix dimension';

    if (nargout == 0)    % Display options.
       disp('pdco default options:')
       disp( defoptions )
    else
       options = defoptions;
    end
    return;
end

Names = ...
[
    'MaxIter    '
    'FeaTol     '
    'OptTol     '
    'StepTol    '
    'StepSame   '
    'x0min      '
    'z0min      '
    'mu0        '
    'Method     '
    'LSQRMaxIter'
    'LSQRatol1  '
    'LSQRatol2  '
    'LSQRconlim '
    'wait       '
    'NOTE       '
];
m     = size (Names,1);
names = lower(Names);

% The remaining clever stuff is from optimset.m.
% We should obtain permission from the MathWorks.

% Combine all leading options structures o1, o2, ... in pdcoSet(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if isstr(arg)                         % arg is an option name
       break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
       if ~isa(arg,'struct')
          error(sprintf(['Expected argument %d to be a ' ...
                'string parameter name ' ...
                'or an options structure\ncreated with pdcoSet.'], i));
       end
       for j = 1:m
          if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
             eval(['val = arg.' Names(j,:) ';']);
          else
             val = [];
          end
          if ~isempty(val)
             eval(['options.' Names(j,:) '= val;']);
          end
       end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end

expectval = 0;                          % start expecting a name, not a value

while i <= nargin
    arg = varargin{i};

    if ~expectval
       if ~isstr(arg)
          error(sprintf(['Expected argument %d to be a ' ...
                         'string parameter name.'], i));
       end

       lowArg = lower(arg);
       j = strmatch(lowArg,names);
       if isempty(j)                       % if no matches
          error(sprintf('Unrecognized parameter name ''%s''.', arg));
       elseif length(j) > 1                % if more than one match
          % Check for any exact matches
          % (in case any names are subsets of others)
          k = strmatch(lowArg,names,'exact');
          if length(k) == 1
             j = k;
          else
             msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
             msg = [msg '(' deblank(Names(j(1),:))];
             for k = j(2:length(j))'
                msg = [msg ', ' deblank(Names(k,:))];
             end
             msg = sprintf('%s).', msg);
             error(msg);
          end
       end
    else
       eval(['options.' Names(j,:) '= arg;']);
    end

    expectval = ~expectval;
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end
    
%
% Copyright (c) 2006. Michael Saunders
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
