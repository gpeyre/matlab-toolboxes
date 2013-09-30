function v2 = prod_vf_sf(v1,s)

% prod_vf_sf - compute the product of a vector field by a scalar field.
%
%   v2 = prod_vf_sf(v1,s)
%
%   The result is the vector field defined by pointwise product.
%
%   Copyright (c) 2004 Gabriel Peyré

v2 = v1 .* repmat(s, [1 1 size(v1,3)]);