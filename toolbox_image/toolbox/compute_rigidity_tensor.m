function H = compute_rigidity_tensor(M,options)

% compute_rigidity_tensor - compute the rigidity
%   tensorial field, which is
%       [(dM/dx)^2    dM/dx*dM/dy]
%       [dM/dx*dM/dy    (dM/dy)^2]
%
%   H = compute_rigidity_tensor(M,options);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'h1')
    options.h1 = 1;
end
h1 = options.h1;
if ~isfield(options, 'h2')
    options.h2 = 1;
end
h2 = options.h2;

if isfield(options, 'h')
    h1 = options.h;
    h2 = options.h;
end

[n,p] = size(M);

grad = compute_grad(M,options);
dx = grad(:,:,1);
dy = grad(:,:,2);

H = zeros(n,p,2,2);
H(:,:,1,1) = dx.^2;
H(:,:,2,2) = dy.^2;
H(:,:,1,2) = dx.*dy;
H(:,:,2,1) = H(:,:,1,2);