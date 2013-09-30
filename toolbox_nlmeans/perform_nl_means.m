function [M1,Wx,Wy,mask_copy] = perform_nl_means(M, options)

% perform_nl_means - denoise an image using non-local Means.
%
%   [M1,Wx,Wy] = perform_nl_means(M,options);
%
%   M is the image to denoise. It can be a color (n,n,s) image
%       where s is the number of chanels.
%   M1 is the denoised image, of same size as M.
%
%   options is a structure that might contains: 
%       k: half size of the neighborhood used to compute windows.
%       T: variance of the gaussian for weights (e.g. 0.05).
%       max_dist: maximum distance to restrict the search.
%       Ma: an image (of same size as M) that provides
%           patches.
%       do_median: set to 0 (default) to perform traditional NLMeans,
%           or set to 1 to perform L1 robust NLMeans 
%           (usefull to deal with salt and pepper noise)
%
%   To avoid manipulating too high dimensional vectors, this code uses
%   a PCA. Set options.ndims to control the dimension of the PCA (e.g. 25).
%
%   Correct the boundary effect by symmetric extension.   
%
%   Copyright (c) 2006 Gabriel Peyre

options.null = 0;
Ma = getoptions(options, 'Ma', M);
if isempty(Ma)
    Ma = M;
end
T = getoptions(options, 'T', .05);
do_median = getoptions(options, 'do_median', 0);
do_patchwise = getoptions(options, 'do_patchwise', 0);
max_dist = getoptions(options, 'max_dist', 10);
mask_process = getoptions(options, 'mask_process', []);
mask_copy = getoptions(options, 'mask_copy', []);
exclude_self = getoptions(options, 'exclude_self', 0);

[m,n,s] = size(M);
[ma,na,sa] = size(Ma);

if isfield(options, 'Vx') && isfield(options, 'Vy')
    Vx = options.Vx;
    Vy = options.Vy;
else
    if na==n && ma==m
        [Vy,Vx] = meshgrid(1:n,1:m);
    else
        Vx = floor( rand(m,n)*(ma-1) ) + 1;
        Vy = floor( rand(m,n)*(na-1) ) + 1;
    end
end

% lift to high dimensional patch space
if not(isfield(options, 'Ha')) || not(isfield(options, 'P')) || not(isfield(options, 'Psi')) 
    [Ha,options.P, options.Psi] = perform_lowdim_embedding(Ma,options);
else
    Ha = options.Ha;
end
if not(isfield(options, 'H'))
    H = perform_lowdim_embedding(M,options);
else
    H = options.H;
end
[M1,Wx,Wy] = perform_nlmeans_mex(Ma,H,Ha,Vx-1,Vy-1,T,max_dist, do_median, do_patchwise, mask_process, mask_copy, exclude_self);
% convert back to matlab notation >0
Wx = Wx + 1;
Wy = Wy + 1;

if not(isempty(mask_process))
    I0 = find(mask_process<0.5); I = [];
    for j=1:s
        I = [I; I0+(j-1)*n*m];
    end
    M1(I) = M(I);
end
if not(isempty(mask_copy))
    mask_copy = 1-double(M1(:,:,1)~=-1);
    M1(M1==-1) = M(M1==-1);
end