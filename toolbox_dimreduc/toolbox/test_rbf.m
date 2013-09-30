% test for linear RBF interpolation
%
%   Copyright (c) 2005 Gabriel Peyré

n = 300;

X = rand(n,2);
f = 0.5 * cos( 5*pi * ( (X(:,1)-0.5).^2 + (X(:,2)-0.5).^2 ) );

p = 30;
x = linspace(0,1,p);
[b,a] = meshgrid(x,x);
Xi = [a(:), b(:)];

fi = perform_rbf_interpolation(Xi,X,f);

M = reshape(fi, p, p);

subplot(1,2,1);
plot_scattered(X,f, 1);

subplot(1,2,2);
imagesc(M);
axis image;