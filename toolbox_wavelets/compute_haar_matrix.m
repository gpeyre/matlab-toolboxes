function M = compute_haar_matrix(w,Jmin, options)

% compute_haar_matrix - compute the haar matrix for a 2D haar transform.
%
%	M = compute_haar_matrix(w,Jmin);
%
% if options.ndims = 2 : 2D transform
%	w is the size of the square (so that M is a transform matrix
%	of size w^2 x w^2).
% if options.ndims = 1 (default) : 1D transform
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;

if isfield(options, 'ndims')
    ndims = options.ndims;
else
    ndims = 1;
end

if exist('perform_haar_transform')
    % use mex file when possible
    haar_callback = @perform_haar_transform;
else
    haar_callback = @perform_haar_transform_slow;
end

if ndims==1
    M = zeros(w,w);
    for i=1:w
        y = dirac(w, i);
        x = feval(haar_callback, y, Jmin, 1);
        M(:,i) = x(:)';
    end
    M = M';
    return;
end

M = zeros(w^2,w^2);
k = 0;
for i=1:w
    for j=1:w
        k = k+1;
        y = dirac([w w], [i,j]);
        x = feval(haar_callback, y, Jmin, -1);
        M(k,:) = x(:)';
    end
end

function y = dirac(s,n)

% dirac - dirac function.
%
%   y = dirac(s,n);
%
%   s is the size of the matrix, and 
%   n is the location of the dirac (default (1,...,1)).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    n = ones(size(s));
end

if length(s)==1
    y = zeros(s,1);
else
    y = zeros(s);
end
y = array_set_val( y, n, 1 );