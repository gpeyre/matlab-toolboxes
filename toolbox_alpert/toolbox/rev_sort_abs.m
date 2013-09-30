function y = rev_sort_abs(x)

% rev_sort_abs - sort by decreasing order of absolute value
%
%   y = rev_sort_abs(x);
%
%   Copyright (c) 2004 Gabriel Peyré

y = sort(abs(x)); y = reverse(y);