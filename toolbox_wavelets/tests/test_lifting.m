options.filter = 'linear'; 
options.filter = '9-7';

%% 1D %%
n = 1024;
x = linspace(0,1,n)';
x = cos(10*pi*x) + (x>.4);

Jmin = 6;

options.ti = 0;
y = perform_lifting_transform(x, Jmin, +1, options);
x1 = perform_lifting_transform(y, Jmin, -1, options);
disp(['Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);

options.ti = 1;
y = perform_lifting_transform(x, Jmin, +1, options);
x1 = perform_lifting_transform(y, Jmin, -1, options);
disp(['Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);


%% test for RWT
options.ti = 1;
options.use_mex = 1;
options.wavelet_vm = 3;
options.wavelet_type = 'daubechies';
y = perform_wavelet_transform(x, Jmin, +1, options);
x1 = perform_wavelet_transform(y, Jmin, -1, options);
disp(['Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);

%% test for wavelab
options.ti = 1;
options.use_mex = 0;
options.wavelet_vm = 3;
options.wavelet_type = 'daubechies';
y = perform_wavelet_transform(x, Jmin, +1, options);
x1 = perform_wavelet_transform(y, Jmin, -1, options);
disp(['Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);


return;

%% 2D %%
n = 256;
Jmin = 3;
x = load_image('lena', n);
x = rescale(x);

options.ti = 0;
y = perform_lifting_transform(x, Jmin, +1, options);
x1 = perform_lifting_transform(y, Jmin, -1, options);
disp(['Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);

options.ti = 1;
y = perform_lifting_transform(x, Jmin, +1, options);
x1 = perform_lifting_transform(y, Jmin, -1, options);
disp(['Error(should be 0)=' num2str(norm(x-x1)/norm(x))]);
