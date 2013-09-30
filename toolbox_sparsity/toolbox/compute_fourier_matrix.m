function W = compute_fourier_matrix(n)

[Y,X] = meshgrid(0:n-1,0:n-1);
W = 1/sqrt(n) * exp( -1i/n *X.*Y );