function M = perform_directional_filtering(M,H,options)

% perform_directional_filtering - perform filtering along some tensor field
%
%   M = perform_directional_filtering(M,H,options);
%
%   H is either a vector field (matrix of size (n,n,2))
%   or a tensor field (matrix of size (n,n,2,2)).
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;
n = size(H,1);

if size(H,4)==1
    % build a tensor field
    e1 = perform_vf_normalization(H);
    e2 = e1(:,:,2:-1:1); e2(:,:,1) = -e2(:,:,1);
    if isfield(options, 'sigma1')
        sigma1 = options.sigma1;
    else
        sigma1 = 8;
    end
    if isfield(options, 'sigma2')
        sigma2 = options.sigma2;
    else
        sigma2 = 3;
    end
    sigma1 = sigma1*ones(size(H,1));
    sigma2 = sigma2*ones(size(H,1));
    H = perform_tensor_recomp(e1,e2,sigma1,sigma2);
elseif size(H,4)~=2
    error('H must be a tensor field or a vector field');
end

if isfield(options, 'n_theta')
    n_theta = options.n_theta;
else
    n_theta = 12;
end
if isfield(options, 'n_sigma')
    n_sigma = options.n_sigma;
else
    n_sigma = 5;
end
if isfield(options, 'n_aniso')
    n_aniso = options.n_aniso;
else
    n_aniso = 4;
end

[e1,e2,sigma1,sigma2] = perform_tensor_decomp(H);
theta = mod( atan2(e1(:,:,2), e1(:,:,1)), pi );
aniso = sigma1./sigma2; 
sigma = sigma1.*sigma2;

if std(aniso(:))<1e-5
    n_aniso = 1;
end
if std(sigma(:))<1e-5
    n_sigma = 1;
end

% compute ranges
s = sort(sigma(:)); a = sort(aniso(:));
k = round(0.05*length(s));
s_list = linspace(s(k), s(end-k), n_sigma);
a_list = linspace(a(k), a(end-k), n_aniso);
t_list = linspace(0,pi,n_theta+1); t_list(end) = [];

% compute the parameters of the filters
[theta_list,sigma_list,aniso_list] = ndgrid(t_list,s_list,a_list);
theta_list = theta_list(:);
aniso_list = aniso_list(:);
sigma1_list = sqrt( sigma_list(:) .* aniso_list );
sigma2_list = sqrt( sigma_list(:) ./ aniso_list );
p = length(theta_list);

% quantize
T = repmat( reshape(theta_list,[1 1 p]), [n n 1] );
S1 = repmat( reshape(sigma1_list,[1 1 p]), [n n 1] );
S2 = repmat( reshape(sigma2_list,[1 1 p]), [n n 1] );
t = repmat( theta, [1 1 p] );
s1 = repmat( sigma1, [1 1 p] );
s2 = repmat( sigma2, [1 1 p] );
E = abs(T-t) + abs(S1-s1) + abs(S2-s2);
[tmp,I] = min(E,[],3);

% do the filtering
F = compute_directional_kernel(sigma1_list,sigma2_list,theta_list,n);
M = perform_adaptive_filtering(M,F,I);