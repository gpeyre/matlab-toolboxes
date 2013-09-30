function M = compute_impulse_noise(M,p,sigma,mu)

% compute_impulse_noise - add impulse noise to an image
%
%	M = compute_impulse_noise(M,p,sigma,mu);
%
%	p is the probability that one pixel is corrupted.
%	sigma is the deviation of the corrupted pixels (uniform noise)
%	mu is the mean
%
%	Copyright (c) 2007 Gabriel Peyre

if nargin<2
	p = 0.05;
end
if nargin<3
	sigma = max(M(:)) - min(M(:)); 
end
if nargin<4
	mu = (max(M(:)) + min(M(:)))/2;
end
n = size(M);
t = randperm(n^2); t = t(1:round(p*n^2));
M = M0;
% Mn(t) = mu + sigma*randn(length(t),1);
M(t) = rand(length(t),1)*sigma + mu;