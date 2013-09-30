function y = array_set_val( x, ind, v )

% array_set_val - set a value in a multidimensional array using a vector for the index
%
%   y = array_set_val( x, ind, v );
%
%   Copyright (c) 2004 Gabriel Peyré

ind = num2cell(ind);
x( ind{:} ) = v;
y = x;