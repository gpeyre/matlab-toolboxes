function T = prod_tf_sf(A,B)

% prod_tf_sf - compute the product of 2 tensor fields
%
%   T = prod_tf_sf(A,B)
%
%   The result is the tensor field defined by pointwise product.
%
%   Copyright (c) 2004 Gabriel Peyré

T = zeros(size(A));

T(:,:,1,1) = A(:,:,1,1)*B(:,:,1,1) + A(:,:,1,2)*B(:,:,2,1);
T(:,:,2,1) = A(:,:,2,1)*B(:,:,1,1) + A(:,:,2,2)*B(:,:,2,1);
T(:,:,1,2) = A(:,:,1,1)*B(:,:,1,2) + A(:,:,1,2)*B(:,:,2,2);
T(:,:,2,2) = A(:,:,2,1)*B(:,:,1,2) + A(:,:,2,2)*B(:,:,2,2);