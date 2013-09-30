function [vlist,A] = perform_fluid_dynamics(v,M,options)

% perform_fluid_dynamics - simulate a viscous fluid
%
%   [vlist,A] = perform_fluid_dynamics(v,M,options);
%
%   v is the initial speed vector
%   M is an initial image to diffuse along the flow
%
%   vlist{i} is the speed at time step i
%   A(:,:,i) is the image diffused at step i
%
%   options.viscosity is the amout of diffusion (# pixels per frame) of the speed
%   options.advspeed is the speed of advection (# pixels per frame) of the speed
%   options.viscosity_texture is the amout of diffusion (# pixels per frame) of the image
%   options.niter is the number of iterations of the simulation
%
%   The solver is the semi-Lagrangian method explained in
%       Jos Stam, "Stable Fluids", In SIGGRAPH 99 Conference Proceedings,
%       Annual Conference Series, August 1999, 121-128.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
sigma = getoptions(options, 'viscosity', 4);
sigma_texture = getoptions(options, 'viscosity_texture', sigma/2);
texture_histo = getoptions(options, 'texture_histo', 'linear');
advspeed = getoptions(options, 'advspeed', 1);
display = getoptions(options, 'display', 1);
niter = getoptions(options, 'niter_fluid', 500);
incompressible = getoptions(options, 'incompressible', 1);

n = size(v,1);

if nargout>2
    A = zeros(n,n,niter);
end

ndisp = max(2,niter/80);

for i=1:niter
    % diffusion
    v = perform_blurring(v, sigma, options);
    % advection   
    w = v;
    for k=1:2
        v(:,:,k) = perform_image_advection(v(:,:,k), advspeed*w, options);
    end
    % make incompressible
    if incompressible>0
        [tmp,vinc] = compute_hodge_decompositon(v,options);
        v = incompressible*vinc + (1-incompressible)*v;
    end
    
    %v1 = perform_vf_normalization(v);
    %v = lambda*v + (1-lambda)*v1;
    
    % advect and blurs M along flow
    if not(isempty(M))
        M = perform_blurring(M, sigma_texture, options);
        if not(isempty(texture_histo))
            M = perform_histogram_equalization(M, texture_histo);
        end
        M = perform_image_advection(M, advspeed*w, options);
    end

    if display && mod(i,ndisp)==1
        clf;
        imageplot(clamp(M)); drawnow;
    end
    
    vlist{i} = v;
    if nargout>1
        A(:,:,i) = clamp(M);
    end
end