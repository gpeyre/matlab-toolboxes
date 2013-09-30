function T = perform_voting_ball(M,sigma,c,p,n)

% perform_voting - perform the voting process for the ball field.
%
%   T = perform_voting_ball(M,sigma,p,n);
%
%   'M' is a 2D function made by the user (1 for features and 0 in the background).
%   
%   Optional:
%   'sigma' and 'c' control the stick field (see 'compute_stick_tf').
%   'n' is the size of the kernel (typically smaller than the size of the image).
%   'p' is the precision for the computation of the ball tensor field.
%
%   See also perform_voting.
%
%   Copyright (c) 2004 Gabriel Peyré




if nargin<2
    sigma = 0.2;
end
if nargin<3
    c = 0.1*sigma;
end
if nargin<4
    p = 32;
end

N = length(M);
n = N/2;
n = floor(n/2)*2+1; % should be odd

% compute the ball voting field
B = compute_ball_tf(N,n,sigma,c,p);

T = zeros( [size(M),2,2] );
for i=1:4   
   % T(:,:,i) = conv2(M,B(:,:,i),'same');              % use inplace acyclic convolution
   T(:,:,i) = perform_convolution(M,B(:,:,i));              % use inplace acyclic convolution
%    T(:,:,i) = perform_convolution(M,kernel(:,:,i));   % use inplace symmetric convolution
end