function L = compute_laplacian_matrix(n,options)

% compute_laplacian_matrix - compute the laplacian matrix
%
%   L = compute_laplacian_matrix(n,options);
%
%   options.ndims set the dimensionality (1 or 2) of the operator.
%   options.bound set the boundary conditions ('per' of 'sym')
%
%   Copyright (c) 2007 Gabriel Peyré

options.null = 0;
if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'per';
end

if isfield(options, 'ndims')
    ndims = options.ndims;
else
    ndims = 2;
end


if ndims==1
    L = -2*speye(n);
    for i=1:n
        L(i, mod(i-2,n)+1) = 1;
        L(i, mod(i,n)  +1) = 1;
    end
    return;
end


% compute the laplacian matrix
L = -4*speye(n^2);

[Y,X] = meshgrid(1:n,1:n);

deltax = [-1 1 0 0];
deltay = [0 0 -1 1];
for i=1:4
    Xi = X + deltax(i);
    Yi = Y + deltay(i);
    if strcmp(bound,'per')
        Xi = mod(Xi-1,n)+1;
        Yi = mod(Yi-1,n)+1;
    else
        I = find(Xi<1); Xi(I) = 2 - Xi(I);
        I = find(Yi<1); Yi(I) = 2 - Yi(I);
        I = find(Xi>n); Xi(I) = 2*n+1 - Xi(I);
        I = find(Yi>n); Yi(I) = 2*n+1 - Yi(I);
    end
    k1 = X + n*(Y-1);
    k2 = Xi + n*(Yi-1);
    k = k1 + n^2*(k2-1);
    L(k(:)) = 1;
end