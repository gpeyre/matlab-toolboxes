function [T,options] = compute_best_threshold(type, M,M0,sigma,options)

% compute_best_threshold - oracle estimation of the threshold
%
% 	T = compute_best_threshold(type, M,M0,sigma,options, G);
%
%	M is the noisy image.
%	M0 is the original image.
%	sigma is the noise level.
%
%	options.G is the geometry (optional).
%
%	Test a range of threshold and find the best.
%	Leave G undefined for wavelet transform.
%	Set G to actual geometry for bandlet transform.
%
%	Copyright (c) 2007 Gabriel Peyre

if nargin<5
    G = []; 
end

if strcmp(type, 'nlmeans')
    %% do denoising
    if isfield(options, 'klist')
        klist = options.klist;
    else
        klist = [3];
    end
    if isfield(options, 'Tlist')
        Tlist = options.Tlist;
    else
        Tlist = [0.02 0.03 0.05 0.07];
        Tlist = linspace(0.02,0.05,10);
    end
    [klist,Tlist] = meshgrid(klist,Tlist);
    err = [];
    for i=1:length(Tlist(:))
        progressbar(i,length(Tlist(:)))
        options.k = klist(i);          % half size for the windows
        options.T = Tlist(i);       % width of the gaussian, relative to max(M(:))  (=1 here)
        [MN1,Wx,Wy] = perform_nl_means(M, options);
        err(end+1) = psnr(M0,MN1);
    end
    [tmp,i] = max(err);
    options.k = klist(i);
    options.T = Tlist(i);
    T = options.T;
    if options.T==min(Tlist(:)) || options.T==max(Tlist(:)) ...
         || options.k==min(klist(:)) || options.k==max(klist(:)) 
        warning('NLmeans Threshold detection: out of bound reached.');
    end
    return;
end


% determine the type
switch lower(type)
    case 'wavelet'
        niter = 25;
        amin = 2; amax = 3.5;
    case 'bandlet'
        if not(isfield(options, 'G'))
            error('You should provide options.G');
        end
        G = options.G;
        niter = 12;
        amin = 2; amax = 3.5;
    case 'pcalet'
        niter = 12;
        amin = 3; amax = 6;
end
Tlist = linspace(amin,amax,niter)*sigma;


err = [];
for i = 1:niter
    T = Tlist(i);
    progressbar(i,niter);

    switch lower(type)
        case 'wavelet'
            % wavelets
            Jmin = 4;
            MW = perform_atrou_transform(M,Jmin,options);
            MWT = perform_thresholding(MW, T);
            M1 = perform_atrou_transform(MWT,Jmin,options);
        case 'pcalet'
            % pcalets
            [MB,options.Pmat] = perform_pcalet_transform( M, options );
            % perform thresholding
            MBT = perform_thresholding(MB, T);
            % reconstruct
            M1 = perform_pcalet_transform( MBT, options );
        case 'bandlet'
            % bandlets
            MB = perform_bandlet_transform( M, G, options );
            % perform thresholding
            MBT = perform_thresholding(MB, T);
            % reconstruct
            M1 = perform_bandlet_transform( MBT, G, options );
    end
    err(i) = psnr(M0,M1);
end

% clf; plot(Tlist/sigma,err); axis tight;

[tmp,I] = max(err);

if I==1 || I==niter
    warning('Out of bound reached.');
end
T = Tlist(I);