function y = callback_atrou(x,dir,options)

% callback_atrou - MCA callback for redundant wavelets.
%
%   y = callback_atrou(x,dir,options);
%
%   Copyright (c) 2007 Gabriel Peyre


if isfield(options, 'Jmin')
    Jmin = options.Jmin;
else
    Jmin = 2;
end

y = perform_atrou_transform(x,Jmin,options);