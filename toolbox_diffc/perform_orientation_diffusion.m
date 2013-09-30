function M = perform_orientation_diffusion(M,t,nsteps,m, bound)

% perform_orientation_diffusion - perform diffusion of orientation data
%
%   M = perform_orientation_diffusion(M,t,nsteps,m);
%
%   M is a scalar matrix that is assumed to contains values modulo m, 
%   (for example m=2pi for data that represent angles, and
%   m=pi for data that represent orientations).
%
%   This function compute the evolution of the heat equation for
%   a total time t using nsteps time steps.
%   But this function treat the evolution as being angle-valued.
%
%   See this publication for more information
%       Perona, "Orientation diffusions"
%       IEEE Transactions on Image Processing, 1998
%   
%
%   Copyright (c) 2005 Gabriel Peyré

if t<0 || nsteps<=0
    return;
end

if nargin<5
    bound = 'sym';
end

% goes into the complex domain
Mx = cos( (2*pi/m) * M );
My = sin( (2*pi/m) * M );

% time step
dt = t/nsteps;

% convolution mask
h = [ 0 dt/4 0; dt/4 1-dt dt/4; 0 dt/4 0 ];

if h(2,2)<0.3
    warning('Instabilities, increase nsteps or decrease t');
end

for i=1:nsteps
    Mx = perform_convolution(Mx,h, bound);
    My = perform_convolution(My,h, bound);
    d = sqrt( Mx.^2 + My.^2 );
    % avoid 0 division
    I = find(d<1e-7); d(I) = 1;
    % project onto unit circle
    Mx = Mx./d;
    My = My./d;
end

M = atan2(My,Mx)*m/(2*pi);