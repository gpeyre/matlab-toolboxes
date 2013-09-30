function [step_types,step_param] = get_lifting_param(type)

% get_lifting_param - get the parameters of the lifting 
%   implementation of a wavelet transform.
%   INTERNAL USE
%
%   [step_types,step_param] = get_lifting_param(type);
%   
%   Copyright (c) 2004 Gabriel Peyré


switch lower(type)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'haar'}
    
    error('Not implemented.');
    % parameter for a 7/9 wavelet transform
    alpha = -1.586134342;
    beta = -0.05298011854;
    gamma = 0.8829110762;
    delta = 0.4435068522;
    zeta = 1.149604398;
    step_types = [0,1,0,1,2];
    step_param = [alpha,beta,gamma,delta,zeta];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'7_9'}
    
    % parameter for a 7/9 wavelet transform
    alpha = -1.586134342;
    beta = -0.05298011854;
    gamma = 0.8829110762;
    delta = 0.4435068522;
    zeta = 1.149604398;
    step_types = [0,1,0,1,2];
    step_param = [alpha,beta,gamma,delta,zeta];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'cubic','4_2'}
    
    % parameter for a cubic (spline) wavelet transform
    % it is the same as the 4/2 wavelet
    step_types = [0,1,0,2];
    step_param = [1/4,1,-3/16,1/2];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'linear','5_3'}
    
    % parameter for a linear (spline) wavelet transform
    % it is the same as the 5/3 wavelet
    step_types = [0,1,2];
    step_param = [-1/2,1/4,sqrt(2)];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'daub4'}
    
    % orthonormal with 3 vanishing moments
    alpha = -0.4122865950;
    beta = -1.5651362796;
    beta1 = 0.352387657;
    gamma = 0.0284590896;
    gamma1 = 0.4921518449;
    delta = -0.3896203900;
    zeta = 1.9182029462;    
    step_types = [6,6,5,5,6,6,5,5,2];
    step_param = [alpha,0,beta1,beta,gamma,gamma1,delta,0,zeta];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
    
    error('Wrong wavelet type.');

end