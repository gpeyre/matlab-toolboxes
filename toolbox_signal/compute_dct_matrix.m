function M = compute_dct_matrix(w, ndims)

% compute_dct_matrix - compute the DCT matrix for a discrete cosine transform.
%
%	M = compute_dct_matrix(w,ndims);
%
%   ndims==1 or 2 for 1D or 2D data.
%	w is the width of the data to transform 
%       (so that M is a transform matrix
%       of size w^ndims x w^ndims).
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    ndims = 2;
end
if length(w)==2
    ndims = 2;
end

if ndims==1
    M = zeros(w,w);
    k = 0;
    for i=1:w
        k = k+1;
        y = dirac(w, i);
        x = idct2(y);
        M(k,:) = x(:)';
    end
else
    if length(w)==1
        w = [w w];
    end
    ww = prod(w);
    M = zeros(ww,ww);
    Jmin = 1;
    k = 0;
    for i=1:w(1)
        for j=1:w(2)
            k = k+1;
            y = dirac(w, [i,j]);
            x = idct2(y);
            M(k,:) = x(:)';
        end
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