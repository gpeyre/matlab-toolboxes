function T = compute_ball_tf(N,n,sigma,c,p)

% compute_ball_tf - compute the voting kernel of a ball tensor.
%
%   T = compute_ball_tf(N,n,sigma,c,p);
%
%   'N' is the total size of the image (supposed to lie in [0,1]x[0,1]).
%
%       The equation for the decay of the field is
%           DF(s,k) = exp(-(s^2+c*k^2)/sigma^2)
%       where 's' is the length of the circle that joins
%       the two points, and 'k' is its curvature.
%
%   Optional:
%   'sigma' control the scale of the voting field.
%   'c' control the ratio between distance and curvature.
%       The higher 'c' is, the narrower the field will be.
%   'n' is the size of the kernel (should be odd).
%
%   'p' is the number of orientation samples needed to compute the field.
%       (the higher the slower the computation).
%
%   See also compute_stick_tf.
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<1
    error('Not enough arguments');
end
if nargin<2
    n = N/2;
end
if nargin<3
    sigma = 0.2;
end
if nargin<4
    c = 0.1*sigma;
end
if nargin<5
    p = 16;
end


n = floor(n/2)*2+1;
nn = (n-1)/2;

T = zeros(n,n,2,2);
for theta = (0:1/p:1-1/p)*2*pi
    v = [cos(theta);sin(theta)];
    B = compute_stick_tf(v,N,n,sigma,c);
    T = T + B;
end
T = T/p;