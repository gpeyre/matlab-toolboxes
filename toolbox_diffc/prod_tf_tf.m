function T = prod_tf_tf(A,B)

% prod_tf_tf - compute the product of 2 tensor fields
%
%   T = prod_tf_tf(A,B)
%
%   The result is the tensor field defined by pointwise product.
%   A tensor field of dimension d is a matrix of dimension (d_1,...,d_k,d,d)
%
%   Copyright (c) 2004 Gabriel Peyre


if size(A,3)==2 && size(A,4)==2 && size(A,5)==1
    % special case for 2D tensor, to speed up a little bit
    T = zeros(size(A));
    T(:,:,1,1) = A(:,:,1,1)*B(:,:,1,1) + A(:,:,1,2)*B(:,:,2,1);
    T(:,:,2,1) = A(:,:,2,1)*B(:,:,1,1) + A(:,:,2,2)*B(:,:,2,1);
    T(:,:,1,2) = A(:,:,1,1)*B(:,:,1,2) + A(:,:,1,2)*B(:,:,2,2);
    T(:,:,2,2) = A(:,:,2,1)*B(:,:,1,2) + A(:,:,2,2)*B(:,:,2,2);
    return;
end

d = size(A);
if d(end-1)~=d(end)
    error('Wrond sizes');
end

A = reshape(A, [prod(d(1:end-2)) d(end-1) d(end)]);

% 3D tensors
T = A*0;
for i=1:d(end)
    for j=1:d(end)
        for k=1:d(end)
            T(:,i,j) = T(:,i,j) + A(:,i,k).*B(:,k,j);
        end
    end
end

T = reshape(T, d);
