% perform_79_transform - biorthogonal 7/9 wavelet transform
%
%   This is the Mex optimized version.
%
%   y = perform_79_transform(x, Jmin, dir, options);
%
%   Perform a 1D in place
%   wavelet transform of 'x' using a 7/9 wavelet.
%   The boundary conditions are handled using symmetric reflexion.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'dir' : direction of the transform (1 : fwd, -1 : bwd)
%
%   Copyright (c) 2005 Gabriel Peyré