function A = build_vandermonde_matrix(x,d)

% build_vandermonde_matrix - build the Vandermonde matrix associdated to the sampling x.
%
%   A = build_vandermonde_matrix(x,d);
%
%   A is a length(x)x(d+1) matrix, with
%   A(i,j) = x(i)^(j-1)
%
%   Copyright (c) 2004 Gabriel Peyré


x = x(:);
n = length(x);

if nargin<2
    d = n-1;
end

A = ones(n,d+1);
for j = 2:d+1
    A(:,j) = x.*A(:,j-1);
end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old code

[J,I] = meshgrid(0:d,1:length(x));
A = x(I).^J;    