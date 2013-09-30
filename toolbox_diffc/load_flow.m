function V = load_flow(name, n, options)

% load_flow - load a predefined flow
%
%   V = load_flow(name, n);
%
%   name can be either 'rot', 'rotinv', 'rotsink', 'rotsrc', 'mixed'
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;
if nargin<2
    n = 100;
end


if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.12;
end

s = 0;
r = 0;
c = [0;0];
p = [0.5;0.5];
        
switch lower(name)
    case 'rot'
        r = 1;
    case 'rotinv'
        r = -1;
    case 'src'
        s = 1;
    case 'sink'
        s = -1;
    case 'rotsink'
        r = 1;
        s = -1;
    case 'rotsrc'
        r = 1;
        s = 1;
    case 'horizontal'
        r = 0;
        s = 0;
        c = [0, 1];
    case 'vertical'
        r = 0;
        s = 0;
        c = [1, 0];
    case 'mixed'
        p = [[0.2;0.2], [0.8;0.8]];
        r = [1,0];  % rotation
        s = [0,1];  % source
        c = [2,0];
    case 'gaussian_rand'
        V = randn(n,n,2);
        h = compute_gaussian_filter(min(60,n/2)*[1 1],sigma,n/2*[1 1]);
        V(:,:,1) = perform_convolution( V(:,:,1), h );
        V(:,:,2) = perform_convolution( V(:,:,2), h );
        if isfield(options, 'normalize')
            V = perform_vf_normalization(V);
        end
        return;
    case 'diffusion_rand'
        nsteps = n;
        m = 2*pi;
        Dt = n*sigma;
        M = rand(n)*m;
        % smoth using intrinsic diffusion
        M = perform_orientation_diffusion(M,Dt,nsteps,m);
        V = zeros(n,n,2);
        V(:,:,1) = cos(M);
        V(:,:,2) = sin(M);
        return;
    case 'phase_rand'
        if isfield(options, 'alpha')
            alpha = options.alpha;
        else
            alpha = 2;
        end
        if isfield(options, 'nbr_turns')
            k = options.nbr_turns
        else
            k = 6;
        end
        M = gen_noisy_image(n,alpha);
        M = 2*pi*k*rescale(M);
        V = zeros(n,n,2);
        V(:,:,1) = cos(M);
        V(:,:,2) = sin(M);
        return;
end

V = compute_flow(p, r, s, c, n);

if isfield(options, 'clamp')
    V = max( min(V,options.clamp), -options.clamp);
end

if isfield(options, 'normalize')
    V = perform_vf_normalization(V);
end