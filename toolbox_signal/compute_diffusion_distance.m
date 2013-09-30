function d = compute_diffusion_distance(x,y)

% compute_diffusion_distance - compute the diffusion distance
%
%   d = compute_diffusion_distance(x,y);
%
%   x and y should be histograms.
%   Use the methode described in 
%
%   Diffusion Distance for Histogram Comparison 
%   H. Ling and K. Okada 
%   IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Vol. I, pp. 246-253, 2006
%
%   Copyright (c) 2006 Gabriel Peyré

z = x(:)-y(:);
% perform cyclic boundary conditions

nsteps = 15;
n = size(x,1);

sigma = 1.4.^(1:nsteps);
sigma = rescale(sigma, 0.001,0.6);
% compute smoothing for varing width
X = zeros(n,nsteps);
for i=1:nsteps
    h = compute_gaussian_filter( round(n/4)*2+1,sigma(i),n);
    X(:,i) = perform_convolution(z,h);
end

% norm for each scale
d = sum(abs(X),1);
% integration though scales
d = sum( d, 2 );