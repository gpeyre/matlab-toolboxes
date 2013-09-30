% test for fluid simulation

path(path, 'toolbox/');

n = 150;

options.bound = 'per';
options.order = 2;  % derivative order for grad/div

name = 'noise';
name = 'lena';

rep = 'results/fluid/';
if not(exist(rep))
    mkdir(rep);
end

dt = .3;
% diffusion per frame
options.viscosity = 2*dt;
% advection per frame
options.advspeed = 1*dt;
% diffusion of the texture
options.viscosity_texture = .5*dt;

options.texture_histo = 'linear';
options.display = 1;
options.niter_fluid = 400;

% initial image
M = load_image(name, 256);
M = crop(M,n);
M = perform_histogram_equalization(M, 'linear');


% initial flow
if strcmp(name, 'lena')
    v = grad(perform_blurring(M,4), options);
    d = sqrt(sum(v.^2,3));
    d = perform_histogram_equalization(d, linspace(.5,1,n^2));
    v = perform_vf_normalization(v);
    v = v .* repmat(d, [1 1 2]);
    v = cat(3, -v(:,:,2), v(:,:,1));
    [tmp,v] = compute_hodge_decompositon(v,options);    
else
    sigma_flow = 30;
    v = perform_blurring(randn(n,n,2), sigma_flow, options);
    v = perform_vf_normalization(v); 
    [tmp,v] = compute_hodge_decompositon(v,options);
end


[vlist,A] = perform_fluid_dynamics(v,M,options);




%% export movie
filename = [rep name '-fluid.avi'];
options.fps = 10;
compute_movie_file(A, filename, options);