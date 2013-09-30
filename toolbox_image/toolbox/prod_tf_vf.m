function vv = prod_tf_vf(t,v)

% prod_vf_vf - compute the product of a tensor field and a vector field.
%
%   vv = prod_tf_vf(t,v)
%
%   The result is the vector field defined by pointwise product.
%
%   WORKS ONLY FOR 2D VECTOR FIELD
%
%   Copyright (c) 2004 Gabriel Peyré

if size(v,3)~=2
    error('Works only for 2D vector fields.');
end

vv = zeros(size(v));
vv(:,:,1) = t(:,:,1,1).*v(:,:,1) + t(:,:,1,2).*v(:,:,2);
vv(:,:,2) = t(:,:,2,1).*v(:,:,1) + t(:,:,2,2).*v(:,:,2);