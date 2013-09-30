function G = compute_periodic_poisson(d, symmetrize)

% compute_periodic_poisson - solve poisson equation 
%
%   G = compute_periodic_poisson(d,symmetrize);
%
%   Solve
%       Delta(G) = d
%   with periodic boundary condition.
%   G has zero mean.
%
%   Set symmetrize=1 (default 0) if the data divergence d is not periodic.
%   This will double the size and consider symmetric extension.
%
%   Copyright (c) 2007 Gabriel Peyre


n = size(d,1);

if nargin==2 && symmetrize
    d = perform_size_doubling(d);
    G = compute_periodic_poisson(d);
    G = G(1:n,1:n);
    return;
end

% solve for laplacian
[Y,X] = meshgrid(0:n-1,0:n-1);
mu = sin(X*pi/n).^2; mu = -4*( mu+mu' );
mu(1) = 1; % avoid division by 0

G = fft2(d) ./ mu; G(1) = 0;
G = real( ifft2( G ) );

%%
function g = perform_size_doubling(g)

g = [g;g(end:-1:1,:,:)];
g = [g,g(:,end:-1:1,:)];