function [M1,MW,MW1] = perform_steerable_matching(M1,M,options)

% perform_steerable_matching - match multiscale histograms
%
% M1 = perform_steerable_matching(M1,M,options);
%
%   M1 is the image to synthesize.
%   M is the exemplar image.
%
%   This function match the histogram of the image and the histogram 
%   of each sub-band of a steerable pyramid.
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
copy_lowpass = getoptions(options, 'copy_lowpass', 1);

if isfield(options, 'color_mode') && strcmp(options.color_mode, 'pca') && ~isfield(options, 'ColorP') && size(M,3)==3
    [tmp,options.ColorP] = change_color_mode(M,+1,options);    
end
rgb_postmatching = getoptions(options, 'rgb_postmatching', 0);

if size(M,3)==3
    options.niter_synthesis = 1;
    if not(isfield(options, 'color_mode'))
        options.color_mode = 'hsv';
    end
    for iter=1:niter_synthesis
        % color images
        M  = change_color_mode(M, +1,options);
        M1 = change_color_mode(M1,+1,options);
        for i=1:size(M,3)
            M1(:,:,i) = perform_steerable_matching(M1(:,:,i),M(:,:,i), options);
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
        [M1(:,:,i),MW,MW1] = perform_steerable_matching(M1(:,:,i),M(:,:,i),options);
    end
    return;
end

if not(isfield(options, 'nb_orientations'))
    options.nb_orientations = 4;
end

nb_orientations = getoptions(options, 'nb_orientations', 1, 1);
Jmin = getoptions(options, 'JminSynthesis', 4);

for iter=1:niter_synthesis
    % spatial equalization
    M1 = perform_histogram_equalization(M1,M);
    % forward transforms
    MW1 = perform_steerable_transform(M1, Jmin, options);
    if iscell(M)
        MW = M;
    else
        MW = perform_steerable_transform(M, Jmin, options);
    end
    % wavelet domain equalization
    MW1 = perform_histogram_equalization(MW1,MW);
    if copy_lowpass
        % copy verbatim the low pass coefficients
        MW1{end} = MW{end};
    end
    % backward transform
    M1 = perform_steerable_transform(MW1, Jmin, options);
    % spatial equalization
    M1 = perform_histogram_equalization(M1,M);
end