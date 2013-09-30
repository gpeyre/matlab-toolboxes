function y = perform_shannon_interpolation(x,eta)

% perform_shannon_interpolation - spectral interpolation
%
%   y = fft_interp(f,factor);
%
%   Copyright (c) 2003 Gabriel Peyré

N = length(x); N0 = (N-1)/2; P =N*eta;
f = fft(x);
f = eta*[f(1:N0+1); zeros(P-N,1); f(N0+2:N)];
y = real( ifft(f) );