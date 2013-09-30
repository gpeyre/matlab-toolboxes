function y = gen_signal_2d(n,alpha)

% gen_signal_2d -  generate a 2D C^\alpha signal of length n x n.
%   gen_signal_2d(n,alpha) generate a 2D signal C^alpha. 
%
%   The signal is scale in [0,1].
%   
%   Copyright (c) 2003 Gabriel Peyré



% new new method

[Y,X] = meshgrid(0:n-1, 0:n-1);

A = X+Y+1;
B = X-Y+n+1;

a = gen_signal(2*n+1, alpha);
b = gen_signal(2*n+1, alpha);
y = a(A).*b(B);
% M = a(1:n)*b(1:n)';

return;


% new method
h = (-n/2+1):(n/2); h(n/2)=1;
[X,Y] = meshgrid(h,h);
h = sqrt(X.^2+Y.^2+1).^(-alpha-1/2);
h = h .* exp( 2i*pi*rand(n,n) );
h = fftshift(h);
y = real( ifft2(h) );

m1 = min(min(y));
m2 = max(max(y));
y = (y-m1)/(m2-m1);

return;

%% old code

y = rand(n,n); 
y = y - mean(mean(y));
for i=1:alpha
    y = cumsum(cumsum(y)')';
    y = y - mean(mean(y));
end
m1 = min(min(y));
m2 = max(max(y));
y = (y-m1)/(m2-m1);