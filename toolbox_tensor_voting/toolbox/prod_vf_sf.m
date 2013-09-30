function v2 = prod_vf_sf(v1,s)

% prod_vf_sf - compute the product of a vector field by a scalar field.
%
%   v2 = prod_vf_sf(v1,s)
%
%   The result is the vector field defined by pointwise product.
%
%   Copyright (c) 2004 Gabriel Peyré

v2 = zeros(size(v1));
v2(:,:,1) = v1(:,:,1).*s;
v2(:,:,2) = v1(:,:,2).*s;