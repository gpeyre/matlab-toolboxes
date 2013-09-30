function imagesc_std(M, epsi)

% imagesc_std - display a rescale version of the image
%
%   imagesc_std(M, epsi);
%
%   The rescaling is computed according to the standard deviation:
%       M -> (M-mean(M))/std(M)
%   Values are clamped to fit in [-epsi,1+epsi].
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    epsi = 1;
end

I = find(isinf(M) | isnan(M));
M(I) = 0;

M = (M-mean(M(:)))/std(M(:))+1;
M = min( max(M,-epsi), 1+epsi);
M(I) = 1;
imagesc(M);