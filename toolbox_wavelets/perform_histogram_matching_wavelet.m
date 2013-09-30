function MW_src = perform_histogram_matching_wavelet(MW_src,MW_tgt, Jmin, options)

% perform_histogram_matching_wavelet - match the histogram of a wavelet transform
%
% Matching of wavelet coefficients only:
%   options.dotransform=0
%   MW_src = perform_histogram_matching_wavelet(MW_src,MW_tgt,Jmin,options);
% Matching of image + wavelet coefficients:
%   options.dotransform=1
%   M_src = perform_histogram_matching_wavelet(M_src,M_tgt,Jmin,options);
%
%   Match the spacial histogram of the image and the histogram of each
%   wavelet sub-band.
%
%   Works also for color images.
%
%   You can set options.use_histomatching=0 if you want to equalize the 
%   kurtosis and skewness of subbands and not their histograms 
%   (works well for natural images).
%   
%   Copyright (c) 2004 Gabriel Peyr?

if nargin>=4
    if ~isstruct(options)
        error('options should be a structure.');
    end
end

options.null = 0;
dotransform = getoptions(options, 'dotransform', 1);
niter_synthesis = getoptions(options, 'niter_synthesis', 1);

if nargin<3
    Jmin = 3;
end

if niter_synthesis>1
    options.niter_synthesis = 1;
    for i=1:niter_synthesis
        MW_src = perform_histogram_matching_wavelet(MW_src,MW_tgt, Jmin, options);
    end
    return;
end

if size(MW_src,3)==3
    options.color_mode = 'pca';
    [MW_src,options.ColorP] = change_color_mode(MW_src,+1,options);
    [MW_tgt,options.ColorP] = change_color_mode(MW_tgt,+1,options);
    for i=1:size(MW_src,3)
        MW_src(:,:,i) = perform_histogram_matching_wavelet(MW_src(:,:,i),MW_tgt(:,:,i),Jmin, options);
    end
    [MW_src,options.ColorP] = change_color_mode(MW_src,-1,options);
    return;
end

if size(MW_src,3)>1
	for i=1:size(MW_src,3)
        MW_src(:,:,i) = perform_histogram_matching_wavelet(MW_src(:,:,i),MW_tgt(:,:,i),Jmin, options);
    end
    return;
end

% for spacial histogram, do not consider absval
options.absval = 0;
options.rows = 0;
options.cols = 0;

if dotransform == 1
    % perform image extension
    n1 = size(MW_src,1);
    n2 = size(MW_tgt,1);
    n = max(n1,n2);
    n = 2^( ceil(log2(n)) );
    MW_src = perform_image_extension(MW_src,n);
    MW_tgt = perform_image_extension(MW_tgt,n);
    options.wavelet_type = 'biorthogonal_swapped';
    options.wavelet_vm = 4;
    % perform spacial matching
    MW_src = perform_matching(MW_src, MW_tgt, options);
    % compute transform
    MW_src = perform_wavelet_transform(MW_src, Jmin, +1, options);
    MW_tgt = perform_wavelet_transform(MW_tgt, Jmin, +1, options);
    % perform coefficients matching
    options.dotransform = 0;
    MW_src = perform_histogram_matching_wavelet(MW_src,MW_tgt,Jmin, options);
    % undo transforms
    MW_src = perform_wavelet_transform(MW_src, Jmin, -1, options);
    MW_tgt = perform_wavelet_transform(MW_tgt, Jmin, -1, options);
    % perform spacial matching
    MW_src = perform_matching(MW_src, MW_tgt, options);
    MW_src = MW_src(1:n1,1:n1);
    return;
end

if size(MW_src,1)~=size(MW_tgt,1)
    error('Wavelets coefficients should be of the same size.');
end

mm = size(MW_src);
Jmax = log2(mm)-1;

for j=Jmax:-1:Jmin
    for q=1:3
        [selx,sely] = compute_quadrant_selection(j,q);
        MW_src(selx,sely) = match_statistics(MW_src(selx,sely), MW_tgt(selx,sely), options);
    end
end

% match low scales
selx = 1:2^Jmin; sely = 1:2^Jmin;
MW_src(selx,sely) = perform_matching(MW_src(selx,sely), MW_tgt(selx,sely), options);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = perform_matching(x,y, options)

% x = perform_histogram_matching(x,y, options);
x = perform_histogram_equalization(x,y, options);

function M = perform_image_extension(M,n)

m = size(M,1);
k = n-m;
while k>size(M,1)
    M = perform_image_extension(M,size(M,1)*2);
    k = k - size(M,1)/2;
end
M = [M; M(end:-1:end-k+1,:)];
M = [M, M(:,end:-1:end-k+1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = match_statistics(x,y, options)


options.null = 0;
if isfield(options, 'use_histomatching')
    use_histomatching = options.use_histomatching;
else
    use_histomatching = 1;
end

if isfield(options, 'nb_bins')
    nb_bins = options.nb_bins;
else
    nb_bins = 100;
end

if use_histomatching
    options.absval = 1;
    x = perform_matching(x, y, options); 
else
    x = perform_kurtosis_equalization(x,y);
end