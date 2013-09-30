function [m,mm,e1,e2] = compute_mutual_information(h)

% compute_mutual_information - compute mutual information of two random variables
%
%   [m,mm,e1,e2] = compute_mutual_information(h).
%
%   h(i,j)=prob(X=i and Y=j)
%
%   m = sum_{i,j} h(i,j) * log2( h(i,j)/(h(i)*h(j))
%   e1 = -sum_i h1(i)*log2(h1(i))
%   e2 = -sum_i h2(i)*log2(h2(i))
%   mm = 2*m/(e1+e2)
%
%   where h1 and h2 are the marginals.
%
%   Copyright (c) Gabriel Peyre

h1 = sum(h,1); h1 = h1/sum(h1);
h2 = sum(h,2); h2 = h2/sum(h2);
hh = h/sum(h(:));
hh(hh<eps)=eps;
h1(h1<eps)=eps;
h2(h2<eps)=eps;
m = hh .* log2( hh ./ (h2*h1) );
m = sum(m(:));

e1 = -sum(h1.*log2(h1));
e2 = -sum(h2.*log2(h2));

mm = 2*m/(e1+e2);
