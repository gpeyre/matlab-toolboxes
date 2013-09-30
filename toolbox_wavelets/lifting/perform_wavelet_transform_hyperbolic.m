function y = perform_wavelet_transform_hyperbolic(x,J,options)

% fwt_hyperbolic - multidimensional hyperbolic (i.e. fully tensorial) wavelet transform
% 
%   y = fwt_hyperbolic(x,J,options);
%
%   perform a hyperbolic transform of the array x, on J scale.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'type' : the type of 1D wavelet transform used (a string), e.g.
%           'peform_79_transform'.
%
%   Copyright (c) 2005 Gabriel Peyré

if J<1
    y = x;
    return; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    options.null = 0;   % force creation of options
end

if isfield(options,'type')   % kind of transform
    type = options.type;
else
    type = 'peform_79_transform';
end

 y = apply_tensorial(x, type, sprintf('%d', J) );

