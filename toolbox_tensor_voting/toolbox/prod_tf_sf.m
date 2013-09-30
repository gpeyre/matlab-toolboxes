function T2 = prod_tf_sf(T1,s)

% prod_tf_sf - compute the product of a tensor field by a scalar field.
%
%   T2 = prod_tf_sf(T1,s)
%
%   The result is the tensor field defined by pointwise product.
%
%   Copyright (c) 2004 Gabriel Peyré

T2 = zeros(size(T1));
for i=1:4
    T2(:,:,i) = T1(:,:,i).*s;
end