function H = compute_hessian(M,options)

% compute_hessian - compute the hessian
%   tensorial field.
%
%   H = compute_hessian(M);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
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
grad = compute_grad(dx,options);
dxx = grad(:,:,1);
dxy = grad(:,:,2);
grad = compute_grad(dy,options);
dyx = grad(:,:,1);
dyy = grad(:,:,2);


H = zeros(n,p,2,2);
H(:,:,1,1) = dxx;
H(:,:,2,2) = dyy;
H(:,:,1,2) = (dxy+dyx)/2;
H(:,:,2,1) = H(:,:,1,2);