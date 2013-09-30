% test for algorithm to minimize TV norm
name = 'lena';
n = 256;
M = load_image(name);
M = rescale(crop(M, n));

options.lambda = .2;

%% compute initial solution
options.niter = 300;
options.method = 'chambolle';
[Mtv,err1,tv1] = perform_tv_denoising(M,options);
options.x0 = Mtv;

%% run the two algorithms
options.niter = 100;
options.method = 'chambolle';
[Mtv,err1,tv1,llist,Err1] = perform_tv_denoising(M,options);
options.method = 'gradient';
[Mtv,err2,tv2,llist,Err2] = perform_tv_denoising(M,options);

%% display residual errors
E = [Err1; Err2]';
plot(log(E)); axis tight;
legend('Chambolle', 'Gradient'); axis tight;