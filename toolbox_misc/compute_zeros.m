function z = compute_zeros(x,y)

% compute_zeros - compute the location of zeros
%
% z = compute_zeros(x,y)
%
%   find the roots z of f(z)=0 where f is sampled at y=f(x)
%
%   Copyright (c) 2007 Gabriel Peyre

I = find( y(1:end-1).*y(2:end) < 0);

d = y(I+1)-y(I); 
J = find(d==0); d(J) = 1;
z = ( x(I).*y(I+1)-x(I+1).*y(I) )./d;
z(J) = (x(I(J))+x(I(J)+1))/2;
