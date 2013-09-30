function x = callback_inpainting(y,dir,options)

% callback_inpainting - callback for inpainting applications
%
%   x = callback_inpainting(y,dir,options);
%
%   You must provide options.mask.
%   The missing pixels (removed) are mask==1.
%
%   Copyright (c) 2008 Gabriel Peyre

mask = getoptions(options, 'mask', [], 1);

x = y;
x(mask==1) = 0;