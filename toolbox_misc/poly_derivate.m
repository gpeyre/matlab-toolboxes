function Q = poly_derivate(P,order)

% poly_derivate - compute the derivative of the polynomial.
% 
%   Q = poly_derivate(P,order);
%
%   P is P(1)+P(2)*X+...+P(n+1)*X^n
%
%   Copyright (c) 2004 Gabriel Peyré

P = reshape(P, prod(size(P)), 1);

if nargin<2
    order=1;
end
if order==0
    return;    
end

Q = P(2:end);
Q = Q.*(1:length(Q))';

if order>1
    Q = poly_derivate(Q,order-1);
end
