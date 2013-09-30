% test for laplacian resolution with periodic boundary conditions

% generate random function
name = 'regular3';
options.bound = 'per';
n = 128;
M = load_image(name, n, options);
M = 1+rescale(M);

[gx,gy] = grad(M,options);
d = div(gx,gy,options);
G = compute_periodic_poisson(d);

% should be zero
norm(G-M+mean(M(:)))