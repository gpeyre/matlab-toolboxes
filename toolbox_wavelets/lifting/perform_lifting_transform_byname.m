function y = perform_lifting_transform_byname(x, Jmin, dir, type, options)

% perform_lifting_transform_byname - wavelet transform via lifting
%
%   y = perform_lifting_transform_byname(x, Jmin, dir, type, options);
%
%   Perform a 1D in place
%   wavelet transform of 'x' using a wavelet specified via string 'type'.
%   The boundary conditions are handled using symmetric reflexion.
%
%   'in place' means that if Lk and Hk are the low-pass
%   and high-pass coefficients at scale k, then the transform
%   is performed this way:
%
%   Original data:  L0 LO LO LO LO LO LO LO
%   1st step:       L1 H1 L1 H1 L1 H1 L1 H1
%   2nd step:       L2 H1 H2 H1 L2 H1 H2 H1
%   3rd step:       L3 H1 H2 H1 H3 H1 H2 H1
%   (continue this until step J<log2(length(x)))
%
%   To come back to Mallat's ordering, aka
%   Original data:  L0 LO LO LO LO LO LO LO
%   1st step:       L1 L1 L1 L1 H1 H1 H1 H1 
%   2nd step:       L2 L2 H2 H2 H1 H1 H1 H1 
%   3rd step:       L3 H3 H2 H2 H1 H1 H1 H1 
%   use the function 'reorder_coefs'.
%
%  'type' is a string containing the type of the transform, and can
%   be either 'haar', '4_2' or 'cubic', '7_9', '5_3' or 'linear', 'daub4'.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb': control verbosity.
%       - 'dir': direction of the transform (1=fwd, -1=bwd)
%   
%   Copyright (c) 2005 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.null = 1;

if nargin<2
    Jmin = 0;
end
if nargin<3
    dir = 1;
end
if nargin<4
    type = '7_9';    
end

if isfield(options,'use_mex')
    use_mex = options.use_mex;
else
    use_mex = 1;
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begining of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[step_type,step_param] = get_lifting_param(type);

if use_mex && exist('perform_lifting_transform')
    y = perform_lifting_transform(x, step_type, step_param, Jmin, dir);
else
    y = perform_lifting_transform_slow(x, step_type, step_param, Jmin, dir, options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
