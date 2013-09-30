function h = compute_laplacian_distribution( name, m, sigma, alpha )

% compute_laplacian_distribution - compute a laplacian distribution
%
% h = compute_laplacian_distribution( name, m, sigma, alpha );
%
%   name can be :
%       'genlaplacian' : exp(-|x|^alpha / sigma)
%       'laplacian' : exp(-|x| / sigma)
%       'genlaplacian0' : only positive entries
%       'laplacian0' : only positive entries
%
%   Copyright (c) Gabriel Peyré 2006

if nargin<3
    alpha = 1;
end

s = length(sigma);
a = length(alpha);
if name(end)=='0'
    % non symmetric
    t = linspace(0,1,m)';
else
    t = linspace(-1,1,m)';
end
sigma = reshape(sigma(:),1,s);
alpha = reshape(alpha(:),1,1,a);
t = repmat( t, [1,s,a] );
sigma = repmat( sigma, [m,1,a] );
alpha = repmat( alpha, [m,s,1] );
h = exp( -abs(t./sigma).^alpha );

% normalize in order to sum to 1
s = sum(h,1); 
I = find(s<eps); s(I) = eps;
s = repmat(s,[m,1,1]);
h = h ./ s;