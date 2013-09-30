function traject = compute_vf_trajectory(points0, vf, options)

% compute_vf_trajectory - compute the trajectories of a vector field
%
%   traject = compute_vf_trajectory(points0, vf, options);
%
%   For each point points0(:,i) in R^2, traject(:,j,i) is the jth point
%   along the trajectory (x(t),y(t)) satisfying
%       d/dt (x,y) = vf(x,y,:)
%
%   This function is very similar to perform_vf_integration.
%   
%   The number points along the trajectory is options.niter
%   The timestep is options.dt
%
%   Copyright (c) 2007 Gabriel Peyr?

options.null = 0;
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 200;
end
if isfield(options, 'dt')
    dt = options.dt;
else
    dt = 0.2;
end

npoints = size(points0,2);
n = size(vf,1);
m = size(vf,2);

traject = zeros(2,niter+1,npoints);
traject(:,1,:) = reshape(points0, [2 1 npoints]);
points = points0;
for i=1:niter
    progressbar(i,niter);
    % evaluate the gradient at the location
    gx = interp2( vf(:,:,1), points(2,:), points(1,:) );
    gy = interp2( vf(:,:,2), points(2,:), points(1,:) );
    points = points + dt * [gx;gy];
    points(1,:) = clamp(points(1,:), 1,n);
    points(2,:) = clamp(points(2,:), 1,m);
    traject(:,i+1,:) = reshape(points, [2 1 npoints]);
end