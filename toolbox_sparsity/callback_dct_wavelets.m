function y = callback_dct_wavelets(x, dir, options)

% callback_dct_wavelets - callback for union of DCT and Wavelets dictionary
%
%   y = callback_dct_wavelets(x, dir, options);
%
%   Copyright (c) 2008 Gabriel Peyre

options.remove_lowfreq = 0;
if not(isfield(options, 'dct_type'))
    options.dct_type = 'redundant';
end
if dir==+1
    y = callback_localdct(x{1}, dir, options) + callback_atrou({x{2:end}}, dir, options);
%    y = y/2;
elseif dir==-1
    yd = callback_localdct(x,dir,options);
    yw = callback_atrou(x,dir,options);
    y = { yd, yw{:}};
else
    error('Pseudo inverse not implemented');
end