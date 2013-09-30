% test for reducing the TV norm
name = 'lena';
n = 256;
M = load_image(name);
M = rescale(crop(M, n));

t = compute_total_variation(M);

tvtgt = t*.1;
niter = 400;

options.tvtgt = tvtgt;
options.niter = niter;
options.nrefresh = 200;
[Mtv,err,tv,lalist] = perform_tv_denoising(M,options);

L = 1/2*err.^2+lalist(end)*tv;
clf;
subplot(2,1,1);
plot(tv-tvtgt, '.-'); title('TV-TV_0'); axis tight;
subplot(2,1,2);
plot(L, '.-'); title('Lagrangian'); axis tight;
