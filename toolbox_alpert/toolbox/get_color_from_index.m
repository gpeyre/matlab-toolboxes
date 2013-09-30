function s = get_color_from_index(i)

% get_color_from_index - turn a number into a color string.
%
%   s = get_color_from_index(i);
%
%   Copyright (c) 2004 Gabriel Peyré

c = {'b' 'g' 'r' 'c' 'm' 'y' 'k'};
s = c{ mod(i-1,7)+1 };
