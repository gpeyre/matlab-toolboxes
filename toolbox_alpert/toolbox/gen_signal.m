function y = gen_signal(n,alpha)

%   gen_signal -  generate a 1D C^\alpha signal of length n.
%
%   y = gen_signal(n,alpha);
%
%   The signal is scale in [0,1].
%   
%   Copyright (c) 2003 Gabriel Peyré



y = randn(n,1); 
fy = fft(y);
fy = fftshift(fy);
% filter with |omega|^{-\alpha}
h = (-n/2+1):(n/2);
h = (abs(h)+1).^(-alpha-0.5);
fy = fy.*h';
fy = fftshift(fy);
y = real( ifft(fy) );

y = (y-min(y))/(max(y)-min(y));

return;

%% old code
y = rand(n,1); 
y = y - mean(y);
for i=1:alpha
    y = cumsum(y);
    y = y - mean(y);
end
y = (y-min(y))/(max(y)-min(y));