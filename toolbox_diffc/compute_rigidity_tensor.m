function H = compute_rigidity_tensor(M,options)

% compute_rigidity_tensor - compute the rigidity tensor field
%
%   H = compute_rigidity_tensor(M,options);
%
%   H is the tensor
%       [(dM/dx)^2    dM/dx*dM/dy]
%       [dM/dx*dM/dy    (dM/dy)^2]
%
%   Copyright (c) 2004 Gabriel Peyre

options.null = 0;


[dx dy] = grad(M,options);

[n,p] = size(M);
H = zeros(n,p,2,2);
H(:,:,1,1) = dx.^2;
H(:,:,2,2) = dy.^2;
H(:,:,1,2) = dx.*dy;
H(:,:,2,1) = H(:,:,1,2);