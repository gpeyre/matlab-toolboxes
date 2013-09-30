function D = compute_synthetic_dictionary(name, w, options)

% compute_synthetic_dictionary - compute a synthetic dictionary
%
%   D = compute_synthetic_dictionary(name, w, options);
%
%   w is the width of the patches.
%   names is 'edges', 'oscillations', 'lines' or 'crossings'.
%
%   Copyright (c) 2007 Gabriel Peyre

options.rescale = 0;
options.manifold_type = name;
switch lower(name)
    case {'edges', 'lines'}
        D = compute_edge_patches(w, options);
    case 'crossings'
        D = compute_crossing_patches(w, options);
    case 'oscillations'
        D = compute_oscillations_patches(w, options);
    otherwise
        error('Unknown type.');
end

p = size(D,3);

% normalize
rescale = getoptions(options, 'rescale', 0);
bell = getoptions(options, 'bell', 'constant');
if rescale
    B = compute_bell_shape(bell, w);
    B = repmat(B,[1,1,p]);
    s = sqrt( sum( sum(B .* (D.^2),1), 2) );
    D = D ./ repmat(s,[w,w,1]);
else
    D = (D+1)/2;
end


% shuffle
D = reshape(D, w^2, p);
D = D(:,randperm(p));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,delta_list,theta_list,freq_list] = compute_oscillations_patches(w, options)

% [M,delta_list,theta_list] = compute_edge_patches(w, options);

options.null = 0;
n_theta = getoptions(options, 'n_theta', 12);
n_freq = getoptions(options, 'n_freq', 6);

t = linspace(0,pi,n_theta+1); t(end) = [];
f = linspace(3,8,n_freq);

[theta_list,freq_list] = meshgrid( t, f );
theta_list = theta_list(:);
freq_list = freq_list(:);
p = length(freq_list);

x = linspace( -(w-1)/2,(w-1)/2, w );
[Y,X] = meshgrid(x,x);
M = []; % zeros(w,w,p*);
% compute images
for k=1:p
    t = theta_list(k);
    f = freq_list(k);
    y = cos(t)*X + sin(t)*Y;
    nbt = round(f*2);
    d = linspace(0,f,nbt+1); d(end)=0;
    for t=1:nbt
        M(:,:,k) = cos( 2*pi* (y+d(t))/f );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,delta_list,theta_list] = compute_edge_patches(w, options)

% [M,delta_list,theta_list] = compute_edge_patches(w, options);

options.null = 0;
n_delta = getoptions(options, 'n_delta', 11);
n_theta = getoptions(options, 'n_theta', 12);
sigma = getoptions(options, 'sigma', .1);
manifold_type = getoptions(options, 'manifold_type', 'edges');

sigma = sigma*w;

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
    switch manifold_type
        case 'edges'
            A = tanh(y/sigma);
        case 'lines'
            A = 2*( 1 - exp( -y.^2 / (2*sigma^2) ) ) - 1;
        case 'oscillations'
            if isfield(options, 'freq')
                freq = options.freq;
            end
            A = cos( y*2*pi/freq );
    end
    M(:,:,k) = A;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, delta1_list,theta1_list, delta2_list,theta2_list] = compute_crossing_patches(w, options)

% [M, delta1_list,theta1_list, delta2_list,theta2_list] = compute_crossing_patches(w, options);

options.null = 0;
n_delta = getoptions(options, 'n_delta', 11);
n_theta = getoptions(options, 'n_theta', 12);
sigma = getoptions(options, 'sigma', .1);

sigma = sigma*w;

dmax = sqrt(2)*w/2+2*sigma;
eta = 2;
t = linspace(0,pi-pi/n_theta,n_theta); t = t(:);
d = linspace(-1,1,n_delta);
d = sign(d).*abs(d).^eta; 
d = d(:)*dmax;

[delta1_list,theta1_list,delta2_list,theta2_list] = ndgrid( d, t, d, t );
theta1_list = theta1_list(:);
delta1_list = delta1_list(:);
theta2_list = theta2_list(:);
delta2_list = delta2_list(:);
p = length(delta1_list);

x = linspace( -(w-1)/2,(w-1)/2, w );
[Y,X] = meshgrid(x,x);

c1 = repmat( reshape(cos(theta1_list),[1 1 p]), [w w 1] );
s1 = repmat( reshape(sin(theta1_list),[1 1 p]), [w w 1] );
c2 = repmat( reshape(cos(theta2_list),[1 1 p]), [w w 1] );
s2 = repmat( reshape(sin(theta2_list),[1 1 p]), [w w 1] );
d1 = repmat( reshape(delta1_list,[1 1 p]), [w w 1] );
d2 = repmat( reshape(delta2_list,[1 1 p]), [w w 1] );
X = repmat( X, [1 1 p] );
Y = repmat( Y, [1 1 p] );

A1 = c1.*X + s1.*Y - d1;
A1 = 2*( 1 - exp( -A1.^2 / (2*sigma^2) ) ) - 1;
A2 = c2.*X + s2.*Y - d2;
A2 = 2*( 1 - exp( -A2.^2 / (2*sigma^2) ) ) - 1;
M = min(A1,A2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = compute_bell_shape(bell, w)

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

