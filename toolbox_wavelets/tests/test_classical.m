% test for the 2D wavelet transform

options.bound = 'per';

g = MakeONFilter('Battle',1);
g = MakeONFilter('Haar',1);
h = mirror_filter(g);

n = 32;
x = (0:1/(n-1):1)';
% x = ones(n,1);
y = x.^2;

Jmin = 1;
f = perform_wavelet_transform_classical(y,Jmin,1,g,h,options);
yy = perform_wavelet_transform_classical(f,Jmin,-1,g,h,options);
