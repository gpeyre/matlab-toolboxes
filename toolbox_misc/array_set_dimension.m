function y = array_set_dimension(x, ind, s, v, range)

% array_set_dimension - Set the value of a vector along a given dimension.
%
%   array_set_dimension(x, ind, s, range) set a vector in matrix 'x' of d
%   dimensions at position 'ind' in dimension 's', i.e. it set
%       x(ind(1),ind(2),...,ind(s-1),:,ind(s),...,ind(d-1)) = v;
%   You can select a subset by using a string 'range', e.g. range='3:6'
%
%   Copyright (c) 2003 Gabriel Peyré

d = length(size(x));

if length(ind)~=(d-1)
   error(['ind should be of dimension d-1. d=', num2str(d), '; ind=',num2str(length(ind)), '.']); 
end
if nargin<4
    error('Not enough argument.')
end
if nargin==4
    range=':';
end

y = x;

% build the corresponding selection string
str_sel = '';
for i=1:(s-1)
    str_sel = [str_sel,int2str(ind(i))];
    str_sel = [str_sel,','];
end;
str_sel = [str_sel,range];
if( s~=d ) str_sel = [str_sel,',']; end;
for i=(s+1):d
    str_sel = [str_sel,int2str(ind(i-1))];
    if( i~=d ) str_sel = [str_sel,',']; end;
end;


if( d==2 && s==2 )
    % special case for 2D array
    str_stock = ['y(',str_sel, ')=v'';'];
else
    str_stock = ['y(',str_sel, ')=v;'];
end
eval(str_stock);