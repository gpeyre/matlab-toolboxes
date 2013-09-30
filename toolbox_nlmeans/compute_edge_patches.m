function [M,delta_list,theta_list] = compute_edge_patches(w, options)

% [M,delta_list,theta_list] = compute_edge_patches(w, options);

options.null = 0;
if isfield(options, 'n_delta')
    n_delta = options.n_delta;
else
    n_delta = 11; 
end
if isfield(options, 'n_theta')
    n_theta = options.n_theta;
else
    n_theta = 12;
end
if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.1;
end
sigma = sigma*w;
if isfield(options, 'rescale')
    rescale = options.rescale;
else
    rescale = 1;
end
% the window to weight the distance
if isfield(options, 'bell')
    bell = options.bell;
else
    bell = 'constant';
end
if isfield(options, 'manifold_type')
    manifold_type = options.manifold_type;
else
    manifold_type = 'edges';
end

switch lower(bell)
    case 'constant'
        B = ones(w);
    case 'sine'
        x = linspace(0,1,w);
        x = (1-cos(2*pi*x))/2;
        B = x'*x;
    otherwise
        error('Unknown bell shape');
end

dmax = sqrt(2)*w/2+2*sigma;
eta = 2;
t = linspace(0,2*pi-2*pi/n_theta,n_theta); t = t(:);
d = linspace(-1,1,n_delta);
d = sign(d).*abs(d).^eta; 
d = d(:)*dmax;


[theta_list,delta_list] = meshgrid( t, d );
theta_list = theta_list(:);
delta_list = delta_list(:);
p = length(delta_list);

x = linspace( -(w-1)/2,(w-1)/2, w );
[Y,X] = meshgrid(x,x);
M = zeros(w,w,p);
% compute images
for k=1:p
    t = theta_list(k);
    d = delta_list(k);
    y = cos(t)*X + sin(t)*Y - d;
    if strcmp(manifold_type, 'edges')
        A = tanh(y/sigma);
    else
        A = 2*( 1 - exp( -y.^2 / (2*sigma^2) ) ) - 1;
    end
    M(:,:,k) = A;
end

% normalize
if rescale
    B = repmat(B,[1,1,p]);
    s = sqrt( sum( sum(B .* (M.^2),1), 2) );
    M = M ./ repmat(s,[w,w,1]);
end