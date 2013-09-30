function x = callback_tomography(y,dir,options)

% callback_tomography - callback for tomography applications
%
%   x = callback_inpainting(y,dir,options);
%
%   You must provide options.mask.
%   mask==1 is the location of KEPT frequencies.
%
%   You must provide options.size which is the size of the signal (signal
%   or image).
%
%   Copyright (c) 2008 Gabriel Peyre

mask = getoptions(options, 'mask', [], 1);

s = getoptions(options, 'size', 1,1);
d = 2; % number of dimensions
if s(1)==1 || s(2)==1
    d = 1;
end

if dir==1
    % measurements
    if d==2
        x = fft2(y) / sqrt( prod(s) );
    else
        x = fft(y) / sqrt( prod(s) );
    end
    x = fftshift(x);
    x = x(mask==1);
else
    % transpose
    x = zeros(s);
    x(mask==1) = y;
    x = fftshift(x);
    if d==2
        x = ifft2(x) * sqrt( prod(s) );
    else
        x = ifft(x) * sqrt( prod(s) );
    end
end