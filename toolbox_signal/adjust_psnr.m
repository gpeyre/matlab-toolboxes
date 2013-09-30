function [sigma,W] = adjust_psnr(W,ptgt,vmax)

% adjust_psnr - adjust noise level to fit psnr
%
%   [sigma,w] = adjust_psnr(W,ptgt,vmax);
%
%   sigma is computed so that psnr(M,M+w*sigma,vmax)=ptgt
%   If vmax=M0 is an image, then sigma is computed so that snr(M0,M0+w*sigma)=ptgt
%
%   Copyright (c) 2008 Gabriel Peyre

if size(vmax)==size(W)
%    20*log10(norm(vmax(:))/norm(sigma*w(:))) = ptgt
    sigma = norm(vmax(:)) / norm(W(:)) * 10^(-ptgt/20);
    return;
end

% 20*log10( vmax / e ) = ptgt
e = sqrt( mean( W(:).^2 ) );
sigma = vmax/e * 10^(-ptgt/20);
W = W*sigma;
