function vertex = perform_normal_displacement(vertex,face,rho,options)

% perform_normal_displacement - perfrom a movement of rho in the normal direction
%
%   vertex = perform_normal_displacement(vertex,face,rho,options);
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

simple_mode = getoptions(options, 'simple_mode', 0);
niter = getoptions(options, 'niter', 5);

if simple_mode
    normals = compute_normal(vertex,face);
    vertex = vertex + rho*normals;
    return;
end

% movement along normal
rho = rho/niter;
for k=1:niter
    normals = compute_normal(vertex,face);
    vertex = vertex + rho*normals;
end