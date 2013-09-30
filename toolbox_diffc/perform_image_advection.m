function M = perform_image_advection(M,v, options)

% perform_image_advection - perform image advection along a flow
%
%   M = perform_image_advection(M,v, options);
%
%   Copyright (c) 2008 Gabriel Peyre

bound = getoptions(options, 'bound', 'per');

n = size(M,1);

[Y,X] = meshgrid(1:n,1:n);
[Y1,X1] = meshgrid(1:n+1,1:n+1);

% extend to avoid boundary problems
M(end+1,:) = M(1,:); M(:,end+1) = M(:,1);

% new coordinates
Xi = X - v(:,:,1);
Yi = Y - v(:,:,2);

if strcmp(bound, 'per')
    Xi = mod(Xi-1,n)+1;
    Yi = mod(Yi-1,n)+1;
elseif strcmp(bound, 'sym')
    Xi(Xi<1) = 2-Xi(Xi<1); Xi(Xi>n) = 2*n-Xi(Xi>n);
    Yi(Yi<1) = 2-Yi(Yi<1); Yi(Yi>n) = 2*n-Yi(Yi>n);
else
    error('Unknown kind of boundaries');
end

% interpolation
M = interp2( Y1,X1,M, Yi,Xi );