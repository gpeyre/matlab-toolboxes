function y  = l2error(x)
% l2error - compute the decreasing of the 
%   L^2 error from orthogonal coefficients.
%
%   y  = l2error(x);
%
%   y[k]^2 = \sum_{i>=k}{ x~[i]^2 }
%   where x~ is x sorted in increasing order.
%
%   Copyright (c) 2003 Gabriel Peyré

x = reshape(x, prod(size(x)), 1);
x = rev_sort_abs(x).^2;
y = sqrt( reverse( cumsum( reverse(x) ) ) );