% test for singular points of vf and tf


rep = 'results/singular-pts/';
if not(exist(rep))
    mkdir(rep);
end

name = 'rot';
name = 'src';
name = 'trissector';
name = 'wedge';
n = 128*2;

% compute vf/tf
switch name
    case {'rot' 'src'}
        U = load_flow(name, n);
        U = perform_vf_normalization(U);
    case {'trissector' 'wedge'}
        pos = [.5 .5];
        T = compute_tensor_field(n, pos,name);
        U = perform_tensor_decomp(T);
end

% compute a cool 2D lic texture
options.isoriented = 0;
M0 = perform_blurring(randn(n),0); % size of the spots
M0 = perform_histogram_equalization( M0, 'linear');
options.histogram = 'linear';
options.dt = 0.8; options.M0 = M0;
options.verb = 1; options.flow_correction = 1;
options.niter_lic = 2;
w = 15;
w = 30;
M = perform_lic(U, w, options);


warning off;
imwrite( rescale(M), [rep name '-texture.png'], 'png' );
warning on;
imageplot(M);