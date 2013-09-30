function [e1,e2,l1,l2] = perform_tensor_eigendecomposition(e1,e2,l1,l2)

% perform_tensor_eigendecomposition - decompose a tensor field
%
% analysis:
%   [e1,e2,l1,l2] = perform_tensor_eigendecomposition(T);
% synthesis:
%   T = perform_tensor_eigendecomposition(e1,e2,l1,l2);
%
%   Copyright (c) 2007 Gabriel Peyre


if nargin==1
    T = e1;
    n = size(T,1);
    if size(T,3)==3 && size(T,4)==1
        T = cat(3, T(:,:,1), T(:,:,3), T(:,:,3), T(:,:,2)); T = reshape(T, [n n 2 2]);
    end
    [e1,e2,l1,l2] = perform_tensor_decomp(T);
elseif nargin==4
    H = perform_tensor_recomp(e1,e2,l1,l2);
    e1 = cat(3, H(:,:,1,1), H(:,:,2,2), H(:,:,2,1) );
else
    error('1 or 4 arguments');
end