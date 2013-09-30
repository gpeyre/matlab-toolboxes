function s = prod_vf_vf(v1,v2)

% prod_vf_vf - compute the dot product of 2 vector field.
%
%   s = prod_vf_vf(v1,v2);
%
%   The result is the scalar field defined by pointwise dot product.
%
%   Copyright (c) 2004 Gabriel Peyré


s = sum(v1.*v2, 3);