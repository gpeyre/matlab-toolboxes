function y = perform_79_transform_slow(x, Jmin, dir, options)

% perform_79_transform_slow - biorthogonal 7/9 wavelet transform
%
%   y = perform_79_transform_slow(x, Jmin, dir, options);
%
%   Perform a 1D in place
%   wavelet transform of 'x' using a 7/9 wavelet.
%   The boundary conditions are handled using symmetric reflexion.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'dir' : direction of the transform (1 : fwd, -1 : bwd)
%
%   Copyright (c) 2003 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    error('Not enough parameters');
end
if nargin<=1
    Jmin = 0;
end
if nargin<=2
    dir = 1;
end
options.null = 0;

if isfield(options, 'use_mex')
    use_mex = options.use_mex;
else
    use_mex = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begining of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = perform_lifting_transform_byname(x, Jmin, dir, '7_9', options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
