function y = perform_lifting_transform_slow(x, step_type, step_param, Jmin, dir, options)

% perform_lifting_transform_slow - perform a lifting transform without relying on Mex file
%
%   This is the Matlab non-optimized version of the function.
%   See 'perform_lifting_transform' for the Mex version.
%
%   y = perform_lifting_transform_slow(x, step_type, step_param, Jmin, dir, options);
%
%   x: the 1D signal to transform in place.
%   step_type: an integer array telling the succession of steps.
%       0 is for 'predict' step : detail channel is computed using
%           predicted value from coarse scale.
%       1 is for 'update' step : coarse channel is computed using
%           value from detail channel to enforce moment conservation.
%       2 is for 'scaling' step : coarse channel is multiplied by some
%           'zeta' factor, while coarse channel is divided by 'zeta'.
%       3 is for 'scaling' on details only.
%       4 is for 'scaling' on coarse only.
%       5 is for 'predict' without feedback
%       6 is for 'update' without feedback
%   step_param : same length as step_types, 1 parameter for each step.
%   Jmin: the coarsest scale. If Jmin=0, then it
%       will perform a full transform (i.e. last coarse value is the mean).
%   dir: +1 for forward transform (default), -1 for backward transform.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'disc' : position of a discontinuity.
%
%   Example of use : 
%       % the parameters of a 7/9 biorthgonal Wavelet
%       alpha = -1.586134342;
%       beta = -0.05298011854;
%       gamma = 0.8829110762;
%       delta = 0.4435068522;
%       zeta = 1.149604398;
%       step_types = [0,1,0,1,2];   
%       step_param = [alpha,beta,gamma,delta,zeta];
%       x = sqrt(1:100);
%       y = perform_lifting_transform_slow(x, step_type, step_param, 3);
%
%   NB: if you don't understand lifting, don't use
%       this function directly, but rather its
%       wrappers like 'perform_lifting_transform_byname'.
%
%   NB: After spliting coarse/detail the signal, i.e. x = [s,d]^T, the types of
%   steps correspond to the matrix multiplications :
%
%   Type 0 :    |1 alpha*(1+z^-1)|
%               |0      1        |
%
%   Type 1 :    |1            0|
%               |beta*(1+z)   1|
%
%   Type 2 :    |zeta   0   |
%               |  0  1/zeta|
%
%   Type 3 :    |1   0  |
%               |0  zeta|
%
%   Type 4 :    |zeta 0|
%               |0    1|
%
%   Type 5 :    |1 alpha+alpha'z^-1)|
%               |0         1        |
%
%   Type 6 :    |1               0|
%               |beta+beta'*z)   1|
%  
%   Copyright (c) 2005 Gabriel Peyré

options.null = 1;

if nargin<2
    [step_type,step_param] = get_lifting_param('7_9');
end
if nargin<4
    Jmin = 0;
end
if nargin<5
    dir = 1;
end

n = length(x);
Jmax = log2(n)-1;
L = Jmax-Jmin+1;

if dir==1
    y = lifting1d(x,L,step_type,step_param,options);
else
    y = lifting1d_inverse(x,L,step_type,step_param,options);
end



function y = lifting1d(x,J,step_types,step_param,options)

% lifting1d - 1D wavelet transform via lifting.
%   INTERNAL USE
%
%   y = lifting1d(x,J,step_types,step_param,offs);
%
%   is the 1D in  place wavelet transform of a vector x.
%
%   Copyright (c) 2003 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    error('Not enough arguments');
end
if length(step_types)~=length(step_param)
    error('Step types and steps parameters must be of the same length.');    
end

if nargin<5
    options.null = 0;   % force creation.    
end

n = length(x);

if n<1
    J=0;
elseif J>ceil(log2(n))
    J=ceil(log2(n));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begining of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=0:J-1
    % select
    sel = 1:2^j:n;  % select the coarse coef from previous step
    xs = x(sel);
    % split
    nn = length(xs);
    s = xs( 1:2:nn );  % coarse
    d = xs( 2:2:nn );  % details
    % the P/U steps
    skip = 0;
    for i=1:length(step_types)
        if skip==1
            skip=0;
        else
            switch step_types(i)
            case 0 % predict
                [s,d] = step_predict(s,d,step_param(i), j, options );
            case 1 % update
                [s,d] = step_update(s,d,step_param(i), j, options );
            case 2 % scale
                [s,d] = step_scale(s,d,step_param(i));
            case 3 % scale on details
                [s,d] = step_scale_d(s,d,step_param(i));
            case 4 % scale on coarse
                [s,d] = step_scale_s(s,d,step_param(i));
            case 5 % predict asymetric
                if (i==length(step_types)) || (step_types(i+1)~=5)
                    error('Problem with asymetric parameters.');
                end
                [s,d] = step_predict_1(s,d,step_param(i),step_param(i+1),j, options);
                skip = 1;    % skip next step
            case 6 % update asymetric
                if (i==length(step_types)) || (step_types(i+1)~=6)
                    error('Problem with asymetric parameters.');
                end
                [s,d] = step_update_1(s,d,step_param(i),step_param(i+1),j, options);
                skip = 1;    % skip next step
            end    
        end
    end
    % reset the coef in correct position
    sel1 = 1:2^(j+1):n;
    x(sel1) = s;
    sel2 = (1+2^j):2^(j+1):n;
    x(sel2) = d;
end

y = x;





function y = lifting1d_inverse(x,J,step_types,step_param,options)

% lifting1d_inverse - 1D inverse wavelet transform via lifting.
%   INTERNAL USE
%
%   lifting1d_inverse(x,J,step_types,step_param,offs);
%
%   is the 1D inplace wavelet transform of a vector x.
%
%   x : the 1D signal to transform in place.
%   J : the coarsest scale. If length(x)==J, then it
%       will perform a full transform (i.e. last coarse value is the mean).
%   step_types : an {0,1,2} array telling the succession of steps.
%       0 is for 'predict' step : detail channel is computed using
%           predicted value from coarse scale.
%       1 is for 'update' step : coarse channel is computed using
%           value from detail channel to enforce moment conservation.
%       2 is for 'scaling' step : coarse channel is multiplied by some
%           'zeta' factor, while coarse channel is divided by 'zeta'.
%       3 is for 'scaling' on details only.
%   step_param : same length as step_types, 1 parameter for each step.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'disc' : position of a discontinuity.
%
%   Example of use : 
%       % the parameters of a 7/9 biorthgonal Wavelet
%       alpha = -1.586134342;
%       beta = -0.05298011854;
%       gamma = 0.8829110762;
%       delta = 0.4435068522;
%       zeta = 1.149604398;
%       step_types = [0,1,0,1,2];   
%       step_param = [alpha,beta,gamma,delta,zeta];
%       x = sqrt(1:100);
%       y = lifting1d(x,3,step_types,step_param);
%  
%   Copyright (c) 2003 Gabriel Peyré

if nargin<4
    error('Not enough arguments');
end
nbr_steps = length(step_types);
if nbr_steps~=length(step_param)
    error('Step types and steps parameters must be of the same length.');    
end

if nargin<5
    options.null = 0;   % force creation.    
end

n = length(x);

if n<1
    J=0;
elseif J>ceil(log2(n))
    J=ceil(log2(n));
end


for j=(J-1):-1:0
    
    % get the coefficients
    sel1 = 1:2^(j+1):n;
    s = x(sel1);
    sel2 = (1+2^j):2^(j+1):n;
    d = x(sel2);
    
    % the P/U steps
    skip = 0;
    for i=nbr_steps:-1:1
        if skip==1
            skip=0;
        else
            switch step_types(i)
            case 0 % predict
                [s,d] = step_predict(s,d,-step_param(i), j, options );
            case 1 % update
                [s,d] = step_update(s,d,-step_param(i), j, options );
            case 2 % scale
                [s,d] = step_scale(s,d,1/step_param(i));
            case 3 % scale on details
                [s,d] = step_scale_d(s,d,1/step_param(i));
            case 4 % scale on coarse
                [s,d] = step_scale_s(s,d,1/step_param(i));
            case 5 % predict asymetric
                if (i==1) || (step_types(i-1)~=5)
                    error('Problem with asymetric parameters.');
                end
                [s,d] = step_predict_1(s,d,-step_param(i-1),-step_param(i),j,options);
                skip = 1;    % skip next step
            case 6 % update asymetric
                if (i==1) || (step_types(i-1)~=6)
                    error('Problem with asymetric parameters.');
                end
                [s,d] = step_update_1(s,d,-step_param(i-1),-step_param(i),j,options);
                skip = 1;    % skip next step
            end
        end    
    end

    % select
    sel = 1:2^j:n;  % select the coarse coef from previous step
    % merge
    nn = length(sel);
    % xs( 1:2:nn ) = s;  % coarse
    % xs( 2:2:nn ) = d;  % details
    x( sel1 ) = s;  % coarse
    x( sel2) = d;  % details
    
end

y = x;

function [s1,d1] = step_predict(s,d,alpha, j, options)

% step_predict - Perform a prediction step in a lifting transform.
%   FOR INTERNAL USE
%
%   [s1,d1] = step_predict(s,d,alpha,j,options);
%   
%   compute the new coarse (s1)
%   and fine (d1) coefficients from the old channels (s and d).
%
%   The parameter is given via 'alpha'.
%   The formula is : d1[k] = d[k] + alpha*(s[k]+s[k+1]).
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%   
%   Copyright (c) 2003 Gabriel Peyré

[s1,d1] = step_predict_1(s,d,alpha,alpha,j, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s1,d1] = step_predict_1(s,d,alpha,alpha1,j, options)

% step_predict - Perform a prediction step in a lifting transform.
%   step_predict(s,d,alpha,j,disc) compute the new coarse (s1)
%   and fine (d1) coefficients from the old channels (s and d).
%   FOR INTERNAL USE
%
%   [s1,d1] = step_predict_1(s,d,alpha,alpha1,j, disc);
%
%   The parameter is given via 'alpha'.
%   The formula is : d1[k] = d[k] + alpha*s[k]+alpha1*s[k+1]).
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%   
%   Copyright (c) 2003 Gabriel Peyré


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    options.null = 0;   % force creation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begining of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls = length(s);
ld = length(d);
d1 = zeros(ld,1);

if ls<ld
    error('s should always be longer than d');    
end

for k=1:ld
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the crossing configuration : 0 = no crossing
    %      d[k-1]  s[k]     d[k]     s[k+1]   d[k+1]
    %   ----*--------+--------*--------+--------*---------
    %                    ^         ^    
    %                   (1)       (2)
    disc_conf = 0;
    p_s_k = 2*(k-1)*2^j + 1;        % position along s channel
    p_d_k = (2*(k-1)+1)*2^j + 1;    % position along d channel
    
    if k+1>ls
        disc_conf = 2;    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computation of s(k)
    if( disc_conf~=1 )
        s_k = s(k);
    else
        % un cross d(k) et s(k)
        s_k = s(k+1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computation of s(k+1)
    if disc_conf~=2 
        % normal case
        s_kk = s(k+1);
    else
        s_kk = s(k);
    end

    d1(k) = d(k) + alpha*s_k + alpha1*s_kk;
end

s1 = s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s1,d1] = step_scale(s,d,zeta)
% perform a scaling step : s1=s*zeta, d1=d/zeta

s1=s*zeta;
d1=d/zeta;


function [s1,d1] = step_scale_d(s,d,zeta)

% step_scale_s - FOR INTERNAL USE
%
% perform a scaling step on detail only : s1=s, d1=d*zeta
%
% perform a scaling step on coarse only : s1=zeta*s, d1=d
%   
%   Copyright (c) 2003 Gabriel Peyré

s1=s;
d1=d*zeta;

function [s1,d1] = step_scale_s(s,d,zeta)

% step_scale_s - FOR INTERNAL USE
%
% [s1,d1] = step_scale_s(s,d,zeta);
%
% perform a scaling step on coarse only : s1=zeta*s, d1=d
%   
%   Copyright (c) 2003 Gabriel Peyré

s1=s*zeta;
d1=d;

function [s1,d1] = step_update(s,d,beta, j, options)

% step_update - Perform an update step in a lifting transform.
%   FOR INTERNAL USE
%
%   step_update(s,d,beta, j, options);
%
%   compute the new coarse (s1)
%   and fine (d1) coefficients from the old channels (s and d).
%
%   The parameter is given via 'alpha'.
%   The formula is : s1[k] = s[k] + beta*(d[k]+d[k-1]).
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%   
%   Copyright (c) 2003 Gabriel Peyré

[s1,d1] = step_update_1(s,d,beta,beta,j, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s1,d1] = step_update_1(s,d,beta,beta1,j, options)

% step_update - Perform an update step in a lifting transform.
%   step_update(s,d,beta, disc,j) compute the new coarse (s1)
%   and fine (d1) coefficients from the old channels (s and d).
%   FOR INTERNAL USE
%
%   [s1,d1] = step_update_1(s,d,beta,beta1,j, disc);
%
%   The parameter is given via 'alpha'.
%   The formula is : s1[k] = s[k] + beta*d[k]+beta1*d[k-1].
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'disc' : position of a discontinuity. The lifting will perform 
%           separatly on d(1:disc) and d(1:disc).
%   
%   Copyright (c) 2003 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    options.null = 0;   % force creation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begining of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls = length(s);
ld = length(d);
s1 = zeros(ls,1);


if ls<ld
    error('s should always be longer than d');    
end

for k=1:ls
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the crossing configuration : 0 = no crossing
    %      s[k-1]  d[k-1]    s[k]     d[k]   s[k+1]
    %   ----*--------+--------*--------+--------*---------
    %                    ^         ^    
    %                   (1)       (2)
    disc_conf = 0;
    p_s_k = 2*(k-1)*2^j + 1;        % position along s channel
    p_d_k = (2*(k-1)+1)*2^j + 1;    % position along d channel
    
    if k==1
        disc_conf = 1;    
    end
    if k>ld       
        disc_conf = 2;    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computation of d(k)
    if( disc_conf~=2 )
        d_k = d(k);
    else
        % un cross s(k) et d(k)
        d_k = d(k-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computation of d(k-1)
    if disc_conf~=1 
        % normal case
        d_kk = d(k-1);
    else
        d_kk = d(k);
    end
    
    s1(k) = s(k) + beta*d_k+beta1*d_kk;
end

d1 = d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
