% perform_lifting_transform - perform a lifting transform without relying on Mex file
%
%   This is the Mex optimized version of the function.
%   See 'perform_lifting_transform_slow' for the Matlab version.
%
%   y = perform_lifting_transform(x, step_type, step_param, Jmin, dir);
%
%   x: the 1D signal to transform in place.
%   step_type: an integer array telling the succession of steps.
%       0 is for 'predict' step : detail channel is computed using
%           predicted value from coarse scale.
%       1 is for 'update' step : coarse channel is computed using
%           value from detail channel to enforce moment conservation.
%       2 is for 'scaling' step : coarse channel is multiplied by some
%           'zeta' factor, while coarse channel is divided by 'zeta'.
%   step_param : same length as step_types, 1 parameter for each step.
%   Jmin: the coarsest scale. If Jmin=0, then it
%       will perform a full transform (i.e. last coarse value is the mean).
%   dir: +1 for forward transform (default), -1 for backward transform.
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
%       y = perform_lifting_transform(x, step_type, step_param, 3);
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
%   Copyright (c) 2005 Gabriel Peyré