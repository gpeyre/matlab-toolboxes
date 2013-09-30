function s = prod_vf_vf(v1,v2)

% prod_vf_vf - compute the dot product of 2 vector field.
%
%   s = prod_vf_vf(v1,v2);
%
%   The result is the scalar field defined by pointwise dot product.
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(v1);

s = zeros(n(1:2));
for i=1:size(v1,3)
    s = s + v1(:,:,i).*v2(:,:,i);
end