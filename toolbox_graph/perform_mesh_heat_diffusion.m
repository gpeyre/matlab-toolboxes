function vertex1 = perform_mesh_heat_diffusion(vertex,face,L,options)

% perform_mesh_heat_diffusion - perform heat diffusion
%
%   vertex1 = perform_mesh_heat_diffusion(vertex,face,L,options);
%
%   L can be a laplacian or a string containing the type of laplacian.
%
%   Each coordinate f=vertex1(i,:) for i=1,2,3 is solution of
%       df/dt=-L*f
%   until a time options.Tmax is reached.
%   The discretization is explicit in time with step size options.dt.
%   options.dt should be small enough (CLF condition) so that the scheme is
%   stable.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'Tmax')
    Tmax = options.Tmax;
else
    Tmax = 3;
end
if isfield(options, 'dt')
    dt = options.dt;
else
    dt = 0.05;
end
if isfield(options, 'verb')
    verb = options.verb;
else
    verb = 1;
end
    

if isstr(L)
    options.symmetrize = 0;
    options.normalize = 1;
    L = compute_mesh_laplacian(vertex,face,L,options);
end

niter = round(Tmax/dt);
vertex1 = vertex;
for i=1:niter
    if verb
        progressbar(i,niter);
    end
    vertex1 = vertex1 - dt*vertex1*L';
end