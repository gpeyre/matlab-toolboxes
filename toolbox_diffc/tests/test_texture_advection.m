% test for texture advection

path(path, '../oriented-patterns-anisotropic-bandlets/');
path(path, '../images/');
path(path, '../images/oriented/');


rep = 'results/advection/';
if not(exist(rep))
    mkdir(rep);
end

name = 'lena';
name = 'noise';

n = 128*2;

M = sum(load_image(name,n),3);
M = rescale(M);

options.bound = 'per';

advspeed = .8;
% perform interpolation between two incompressible fluids
sigma_flow = 100;
v = perform_blurring(randn(n,n,2), sigma_flow, options);
[tmp,v] = compute_hodge_decompositon(v,options);
v = perform_vf_normalization(v);

niter = 100;
A = zeros(n,n,2);
for i=1:niter
    A(:,:,i) = M;
    M = perform_image_advection(M, advspeed*v, options);
    imageplot(M); drawnow;
end



%% export movie
filename = [rep name '-advection.avi'];
fps = 10;
compute_movie_file(A, filename, fps);
