function T = compute_gradient_tensor(M,h,options)

% compute_gradient_tensor - compute the structure tensor
%
%   T = compute_gradient_tensor(M,h,options);
%
%   The main eigenvector of T is aligned (approx.) with the gradient
%   direction.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if iscell(M)
    for i=1:length(M)
        T{i} = compute_gradient_tensor(M{i},h,options);
    end
    return
end
if size(M,3)>1
    u = M;
else
    u = grad(M, options);
    % u = compute_grad(M);
end
T = cat(3, u(:,:,1).^2,  u(:,:,2).^2, u(:,:,1).*u(:,:,2) );
T = perform_convolution(T,h, options);