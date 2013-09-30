function c = gen_rand_curve(n, alpha, options)

% gen_rand_curve - generate a curve with C^alpha smoothness
%
%   c = gen_rand_curve(n, alpha, options);
%
%   Copyright (c) 2005 Gabriel Peyré

n = round(n/2)*2;

options.null = 0;
if isfield(options, 'rmax')
    rmax = options.rmax;
else
    rmax = 0.4;
end

% amplitud deacy
h = (-n/2+1):(n/2);
h = (abs(h)+1).^(-alpha-0.5);
h = rescale( h, 0,rmax );
% add random phase
y = h .* exp( rand(1,n)*2i*pi );

y = fftshift(y);
y(1) = 0;   % zero mean
% around a circle
y(2) = 1;
y(n) = 0;

y = ifft(n*y);

% close the curve
y(end+1) = y(1);

c = [real(y); imag(y)];