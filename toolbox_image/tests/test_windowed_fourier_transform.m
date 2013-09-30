% test for local fourier transform

n = 128;
name = 'barb';
M = crop(load_image(name),n);

w = 15; q = 7;
options.window_type = 'sin';
MF = perform_windowed_fourier_transform(M,w,q,n, options);
M1 = perform_windowed_fourier_transform(MF,w,q,n, options);

e = sum(M(:).^2);
ef = sum(abs(MF(:)).^2);
disp( ['Energy conservation error (should be 0) : ' num2str( (e-ef)/e*100 ) '%' ]);