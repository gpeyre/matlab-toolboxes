function y = poly_val(P,x)

% poly_val - evaluate a polynomial
%
%   y = poly_val(P,x);
%
%   P is P(1)+P(2)*X+...+P(n+1)*X^n
%
%   Copyright (c) 2004 Gabriel Peyré

if size(P,1)==1 && size(P,2)~=1
    P = P';
end

if size(x,1)==1 && size(x,2)~=1
    x = x';
end


y = zeros( size(x) );
for i=1:size(x,2)
    y(:,i) = polyval(reverse(P),x(:,i));
end