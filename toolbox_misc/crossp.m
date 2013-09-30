function z = crossp(x,y)

% crossp - compute cross product
%
%   z = crossp(x,y);
%
% x and y are (m,3) dimensional
%
%   Copyright (c) 2007 Gabriel Peyre

z = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);