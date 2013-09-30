function M = compute_tensorial_operator(T,n,p)

% compute_tensorial_operator - compute a 2D convolution operator that act in a tensorial way
%
%   M is a n²xn² matrix that represent an operator
%       M : (R^n x R^p)  ---> (R^n x R^p)
%
%   Example for the discretized Laplacian :
%
%   T = [0 -1 0; -1 4 -1; 0 -1 0];
%   M = compute_tensorial_operator(T,5);
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<3
    p = n;
end


k = size(T,1);
if size(T,2)~=k || mod(k,2)~=1
    errro('T should be a square matrix of odd size');
end
k = (k-1)/2;

[dY,dX] = meshgrid(-k:k,-k:k);
[Y,X] = meshgrid(1:p,1:n);

M = sparse(n*p,n*p);

for s=1:(2*k+1)^2
   dx = dX(s);
   dy = dY(s);
   v = T(s);
   
   
   Z = X(:) + n*( Y(:)-1 );
   Z1 = X(:)+dx + n*( Y(:)+dy-1 );
   
   % remove out of bound points
   I = find( X(:)+dx>n | X(:)+dx<1 | Y(:)+dy>n | Y(:)+dy<1 );
   Z1(I) = [];  Z(I) = [];

   M = matrix_sampling_set(M, v, [Z';Z1']);
   
   if 1
   I = find( (X(:)+dx>n | X(:)+dx<1) & Y(:)+dy<=n & Y(:)+dy>0 );
   Z = X(I) + n*( Y(I)-1 );
   Z1 = X(I)-dx + n*( Y(I)+dy-1 );
   M = matrix_sampling_add(M, v, [Z';Z1']);
   
   I = find( (Y(:)+dy>n | Y(:)+dy<1) & X(:)+dx<=n & X(:)+dx>0 );
   Z = X(I) + n*( Y(I)-1 );
   Z1 = X(I)+dx + n*( Y(I)-dy-1 );
   M = matrix_sampling_add(M, v, [Z';Z1']);
   end
   
end