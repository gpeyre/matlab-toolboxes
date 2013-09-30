function y = compute_orthogonal_projection( x, A, C )

% compute_orthogonal_projection - computer orthogonal projection
%
%   y = compute_orthogonal_projection( x, A, C );
%
%   y is the pojection of x onto the intersection of the sets
%       { y \ <x-cj,aj>=0 }
%   where the cj and the aj are the columns of C and A.
%
%   Copyright (c) 2007 Gabriel Peyre

% compute solution as
%   y = x + A*lambda

D = (A'*A);
if abs(det(D))<1e-9
    D = D+rand(size(D))*1e-9;
end
lambda = D \ ( -A'*x + diag(A'*C) );
y = x + A*lambda;

