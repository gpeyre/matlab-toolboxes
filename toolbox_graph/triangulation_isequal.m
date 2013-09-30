function b = triangulation_isequal(face,faced)

% triangulation_isequal - true if two triangulations are equals
%
%   b = triangulation_isequal(face,faced)
%
%   Copyright (c) 2008 Gabriel Peyre

I = setdiff(sort(face)',sort(faced)', 'rows');
b = isempty(I);