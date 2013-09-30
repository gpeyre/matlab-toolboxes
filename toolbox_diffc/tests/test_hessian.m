% test for hessian tensor

n = 256;
name = 'barb';
M = load_image(name);
M = rescale(crop(M,n));

options.order = 2;
T = compute_hessian(M, options);

[e1,e2,l1,l2] = perform_tensor_decomp(T);
l1 = abs(l1); l2 = abs(l2);
T = perform_tensor_recomp(e1,e2,l2,l1);

sigma = 5;
T = perform_blurring(T,sigma);

options.sub = 8;
clf;
plot_tensor_field(T, M, options);
figure;
clf;
plot_tensor_field_nb(T, M, options);