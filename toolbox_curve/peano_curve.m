function c = peano_curve(n)

% hilbert_curve - create a peano curve at level n of recursion
%
%   c = peano_curve(n);
%
%   Copyright (c) 2004 Gabriel Peyré

N = 100000;
t = 0:1/(N-1):1;

c = zeros(2,N);
for i=1:n
    c(1,:) = c(1,:) + f( 10^(i-1)*t )/2^i;
    c(2,:) = c(2,:) + g( 10^(i-1)*t )/2^i;
end

function y = f(t)
t = mod(t,1);

y = (t>=0.5 & t<0.8)  + (t>=0.4 & t<0.5).*(t-0.4)*10 + (t>=0.8 & t<0.9).*(0.9-t)*10;
% y = (t>=0.5 & t<0.8)  + (t>=0.4 & t<0.5).*(t-0.4)*10 + (t>=0.8 & t<1).*(1-t)*5;

function y = g(t)
t = mod(t,1);

y = (t>=0.3 & t<0.4) + (t>=0.7 & t<0.8) + (t>=0.2 & t<0.3).*(t-0.2)*10 + (t>=0.6 & t<0.7).*(t-0.6)*10 + ...
    (t>=0.4 & t<0.5).*(0.5-t)*10 + (t>=0.8 & t<0.9).*(0.9-t)*10;
% y = (t>=0.3 & t<0.4) + (t>=0.7 & t<0.8) + (t>=0.2 & t<0.3).*(t-0.2)*10 + (t>=0.6 & t<0.7).*(t-0.6)*10 + ...
%    (t>=0.4 & t<0.5).*(0.5-t)*10 + (t>=0.8 & t<1).*(1-t)*5;