function [M1,MW,MW1] = perform_wavelet_matching(M1,M,options)

% perform_wavelet_matching - match multiscale histograms
%
% M1 = perform_wavelet_matching(M1,M,options);
%
%   M1 is the image to synthesize.
%   M is the exemplar image.
%
%   This function match the histogram of the image and the histogram 
%   of each sub-band of a wavelet-like pyramid.
%
%   To do texture synthesis, one should apply several time this function.
%   You can do it by setting the value of options.niter_synthesis.
%   This leads to the synthesis as described in 
%
%       Pyramid-Based Texture Analysis/Synthesis
%       D. Heeger, J. Bergen,
%       Siggraph 1995
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
niter_synthesis = getoptions(options, 'niter_synthesis', 1);

if not(isfield(options, 'color_mode'))
    options.color_mode = 'pca';
end
if isfield(options, 'color_mode') && strcmp(options.color_mode, 'pca') && ~isfield(options, 'ColorP') && size(M,3)==3
    [tmp,options.ColorP] = change_color_mode(M,+1,options);    
end
rgb_postmatching = getoptions(options, 'rgb_postmatching', 0);

if size(M,3)==3
    options.niter_synthesis = 1;
    for iter=1:niter_synthesis
        % color images
        M  = change_color_mode(M, +1,options);
        M1 = change_color_mode(M1,+1,options);
        for i=1:size(M,3)
            M1(:,:,i) = perform_wavelet_matching(M1(:,:,i),M(:,:,i), options);
        end
        M  = change_color_mode(M, -1,options);
        M1 = change_color_mode(M1,-1,options);
        if rgb_postmatching
            for i=1:size(M,3)
                M1(:,:,i) = perform_histogram_equalization(M1(:,:,i),M(:,:,i));
            end
        end
    end
    return;
end

if size(M,3)>1
    for i=1:size(M,3)
        [M1(:,:,i),MW,MW1] = perform_wavelet_matching(M1(:,:,i),M(:,:,i),options);
    end
    return;
end

n = size(M,1);
n1 = size(M1,1);

synthesis_method = getoptions(options, 'synthesis_method', 'steerable');

m = 2^( ceil(log2(n)) );
m1 = 2^( ceil(log2(n1)) );
M = perform_image_extension(M,m);
M1 = perform_image_extension(M1,m1);

% precompute input
MW = my_transform(M, +1, options);

for iter=1:niter_synthesis
    % spatial equalization
    M1 = my_equalization(M1,M);
    % forward transforms
    MW1 = my_transform(M1, +1, options);
    % wavelet domain equalization
    MW1 = my_equalization(MW1,MW);
    % backward transform
    M1 = my_transform(MW1, -1, options);
    % spatial equalization
    M1 = my_equalization(M1,M);
end

M1 = M1(1:n1,1:n1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = my_equalization(M,M0)

options.absval = 0;
options.rows = 0;
options.cols = 0;
options.dim3 = 1;

if iscell(M)
    for i=1:min(length(M),length(M0))
        M{i} = my_equalization(M{i},M0{i});        
    end
    return;
end
if size(M,3)>1
    for i=1:min(size(M,3),size(M0,3))
        M(:,:,i) = my_equalization(M(:,:,i),M0(:,:,i));        
    end
    return;
end
M = perform_histogram_equalization(M,M0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = my_transform(M, dir, options)

n = size(M,1);
deltaJ = 5;
Jmax = log2(n)-1;
Jmin = max(Jmax-deltaJ+1,3);
    
synthesis_method = getoptions(options, 'synthesis_method', 'steerable');
% steerable options
if not(isfield(options, 'nb_orientations'))
    options.nb_orientations = 4;
end
% wave ortho options
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;

switch synthesis_method
    case 'steerable'
        M = perform_steerable_transform(M, Jmin, options);
    case 'wavelets-ortho'
        if dir==-1
            M = convert_wavelets2list(M, Jmin);
        end
        M = perform_wavelet_transform(M, Jmin, dir, options);
        if dir==1
            M = convert_wavelets2list(M, Jmin);            
        end    
    case 'quincunx-ti'
        M = perform_quicunx_wavelet_transform_ti(M,Jmin,options);
    case 'wavelets-ti'
        options.wavelet_type = 'biorthogonal';
        options.wavelet_vm = 3;
        M = perform_atrou_transform(M,Jmin,options);
    otherwise 
        error('Unknown transform.');
end