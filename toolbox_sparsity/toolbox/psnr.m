function p = psnr(x,y, vmax)

% psnr - compute the Peack Signal to Noise Ratio, defined by :
%       PSNR(x,y) = 10*log10( max(max(x),max(y))^2 / |x-y|^2 ).
%
%   p = psnr(x,y);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    m1 = max( abs(x(:)) );
    m2 = max( abs(y(:)) );
    vmax = max(m1,m2);
end

d = mean( (x(:)-y(:)).^2 );

p = 10*log10( vmax^2/d );