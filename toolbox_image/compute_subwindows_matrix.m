function X = compute_subwindows_matrix(A,k,mask)

% compute_subwindows_matrix - compute a matrix containing each sub-windows.
%   
%   X = compute_subwindows_matrix(A,k,mask);
%
%   If A is of size [n,p], then X is of size [(2*k+1)^2,n*p].
%   Each row contains a sub-windows of A of size (2*k+1)x(2*k+1).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    k = 3;
end
if nargin<3
    mask=ones(size(A));
end

mmask = zeros(2*k+size(mask));
mmask(k+1:end-k,k+1:end-k) = mask;
A = symmetric_extension(A,k);

I = find(mmask>0);
p = length(I);
X = zeros((2*k+1)^2,p);
for ii = 1:p
    i = I(ii);
    [i1,i2] = ind2sub(size(A),i);
    Ai = A(i1-k:i1+k, i2-k:i2+k);
    X(:,ii) = Ai(:);
end