function y = gen_brownian_motion(n,d,sigma)

% gen_brownian_motion_1d - generate a 1D brownian motion
%
%   y = gen_brownian_motion(n,d,sigma);
%
%   n is the length, d the dimension, sigma the variance of time increment
%   (should be dt^0.5).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<1
    n = 1000;
end
if nargin<2
    d = 1;
end
if nargin<3
    sigma=1;
end

% for random walk
z= 2.*(rand(d,n)<=0.5)-1; 
% for brownian
z = sigma.*randn(d,n);

y = [zeros(1,d); cumsum(z')];
