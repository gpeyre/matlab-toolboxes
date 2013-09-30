function H = compute_hessian(M,options)

% compute_hessian - compute the hessian tensor field.
%
%   H = compute_hessian(M, options);
%
%   Copyright (c) 2004 Gabriel Peyre

options.null = 0;


if size(M,3)==1
    % 2D field
    [n,p] = size(M);

    [dx,dy] = grad(M,options);
    [dxx,dxy] = grad(dx,options);
    [dyx,dyy] = grad(dy,options);

    H = zeros(n,p,2,2);
    H(:,:,1,1) = dxx;
    H(:,:,2,2) = dyy;
    H(:,:,1,2) = (dxy+dyx)/2;
    H(:,:,2,1) = H(:,:,1,2);
else
    % 3D field
    [n,p,q] = size(M);

    G = grad(M,options);
    H = zeros(n,p,q,3,3);
    for i=1:3
        [H(:,:,:,i,1),H(:,:,:,i,2),H(:,:,:,i,3)] = grad(G(:,:,:,i),options);
    end

end