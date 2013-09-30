% test for the circle haar transform
%   
%   Copyright (c) 2006 Gabriel Peyré

n = 128;
Jmin = 0;
m = pi;
options.normalize = 1;

k = 4;
G = load_image('lena',n);
x = linspace(-1,1,n);
[Y,X] = meshgrid( x,x );
G = X.*Y; % X.*Y.^2 + X.^2;
G = X.*Y;

G = mod( rescale(G)*k*m, m);
% G = round( rand(n,1)*256 )*0+1;

% FWD
G1 = perform_circle_haar_transform(G, Jmin, 1, m, options);
G2 = perform_haar_transform(G, Jmin, 1);
G2r = reorder_coefs(G2,0,1);

% BWD
GG1 = perform_circle_haar_transform(G1, Jmin, -1, m, options);
GG2 = perform_haar_transform(G2, Jmin, -1);

% should be 0 for correct reconstruction
norme( cos(2*GG1) - cos(2*G) ) + norme( sin(2*GG1) - sin(2*G) )
norme( cos(2*GG2) - cos(2*G) ) + norme( sin(2*GG2) - sin(2*G) )

clf;
subplot(2,2,1);
imagesc(G);
title('Original signal');
axis image; axis off;
subplot(2,2,3);
plot_wavelet( G1, Jmin );
title('Circle haar transform');
subplot(2,2,4);
plot_wavelet( G2r, Jmin );
title('Haar transform');
colormap gray(256);

err_circle = [];
nbr_coefs_circle = [];
err_haar = [];
nbr_coefs_haar = [];
T_list = [1:10:1000]*m/256;
for T=T_list
    G1T = G1 .* (abs(G1)>=T);
    G2T = G2 .* (abs(G2)>=T);
    nbr_coefs_circle = [nbr_coefs_circle, sum(abs(G1(:))>=T)];
    nbr_coefs_haar = [nbr_coefs_haar, sum(abs(G2(:))>=T)];
    GG1T = perform_circle_haar_transform(G1T, Jmin, -1, m, options);
    GG2T = perform_haar_transform(G2T, Jmin, -1);
    err_circle = [err_circle, norme( cos(2*GG1T) - cos(2*G) ) + norme( sin(2*GG1T) - sin(2*G) )];
    err_haar =   [err_haar,   norme( cos(2*GG2T) - cos(2*G) ) + norme( sin(2*GG2T) - sin(2*G) )];
end

subplot(2,2,2);
loglog(nbr_coefs_circle,err_circle,'.-',nbr_coefs_haar,err_haar,'.-');
axis tight;
legend('circlehaar', 'haar');
% L^1 norms
n1 = sum(abs(G1(:)))/n^2;
n2 = sum(abs(G2(:)))/n^2;
title( sprintf('L^1(circlehaar)=%.2f  L^1(haar)=%.2f', n1,n2) );

return;


if ~exist(G)
    n = 128;
    [Y,X] = meshgrid(1:n,1:n);
    G = (X+Y);
    G = rescale(G)*pi + (rand(n)>0.5)*pi;
    
end