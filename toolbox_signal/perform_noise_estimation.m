function sigma = perform_noise_estimation(x)

% perform_noise_estimation - estimate additive noise level
%
%   sigma = perform_noise_estimation(x);
%
%   Estimate the std of a gaussian additive noise perturbing
%   an image or a 1D signal.
%   Use the median of high pass wavelet coefficients.
%
%   Copyright (c) 2007 Gabriel Peyre

if size(x,3)>1
    % color images
    for i=1:size(x,3)
        sigma(i) = perform_noise_estimation(x(:,:,i));
    end
    sigma = mean(sigma);
    return;
end

if size(x,1)==1
    x = x';
end

n = size(x);
% compute high pass coefficient of a Haar transform
xd = x;
xd = xd(1:2:end,:) - xd(2:2:end,:); xd = xd' / sqrt(2);
if n(2)>1
    % 2D signals
    xd = xd(1:2:end,:) - xd(2:2:end,:); xd = xd'/ sqrt(2);
end

% mad estimator
sigma = mad(xd(:),1)/0.6745;