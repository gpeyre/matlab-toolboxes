function w = prod_tf_vf(T,v)

% prod_vf_vf - compute the product of a tensor field and a vector field.
%
%   w = prod_tf_vf(T,v)
%
%   The result is the vector field defined by pointwise product of tensor time vector.
%   A tensor field of dimension d is a matrix of dimension (d_1,...,d_k,d,d)
%   A vector field of dimension d is a matrix of dimension (d_1,...,d_k,d)
%
%   Copyright (c) 2004 Gabriel Peyré

if size(v,3)==2 && size(v,4)==1
    % special case for 2D tensor, to speed up a little bit
    w = zeros(size(v));
    w(:,:,1) = T(:,:,1,1).*v(:,:,1) + T(:,:,1,2).*v(:,:,2);
    w(:,:,2) = T(:,:,2,1).*v(:,:,1) + T(:,:,2,2).*v(:,:,2);
    return;
end

d = size(T);
if d(end-1)~=d(end)
    error('Wrond sizes');
end

T = reshape(T, [prod(d(1:end-2)) d(end-1) d(end)]);
v = reshape(v, [prod(d(1:end-2)) d(end-1)]);


% 3D tensors
w = v*0;
for i=1:3
    for j=1:3
        w(:,i) = w(:,i) + T(:,i,j).*v(:,j);
    end
end
w = reshape(w, d(1:end-1));