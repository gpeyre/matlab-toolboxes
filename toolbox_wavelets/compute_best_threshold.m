function [T,M1,pval,Tlist,err] = compute_best_threshold(type, M, M0, sigma, options)

% compute_best_threshold - oracle estimation of the threshold for denoising
%
% 	[T,M1,pval,Tlist,err] = compute_best_threshold(type, M, M0, sigma, options);
%
%	Test a range of thresholds and find the best given a target clean image.
%
%	M is the noisy image.
%	M0 is the original image.
%	sigma is the noise level.
%
%	options.thresh is either 'soft' or 'hard' or 'block'
%
%	T is the best threshold.
%	pval is the corresponding denoising SNR.
%	M1 is the denoised image.
%   Tlist is the list of tested thresholds.
%
%   type is either 'wavti' (translation invariant wavelets) or 
%   'wavortho' (orthogonal wavelets) or 'curvelet'
%
%   Set options.thresh_low=1 to also threshold the low wavelet frenquencies
%   (this is usually not very good for soft thresholding).
%
%   options.nthresh gives the number of tested thresholds.
%
%	Copyright (c) 2008 Gabriel Peyre

options.null = 0;

thresh = getoptions(options, 'thresh', 'hard');
verb = getoptions(options, 'verb', 1);
thresh_low = getoptions(options, 'thresh_low', 0);
% for bandlets
G = getoptions(options, 'geometry', []);

% options for wavelets
Jmin = getoptions(options, 'Jmin', 4);
if not(isfield(options, 'wavelet_type'))
    options.wavelet_type = 'biorthogonal_swapped';
end
if not(isfield(options, 'wavelet_vm'))
    options.wavelet_vm = 4;
end

% number of tested thresholds
niter = getoptions(options, 'nthresh', 15);

% determine the type
switch lower(type)
    case {'wavelet' 'waveti' 'wavti'}
        amin = 2.3; amax = 3.2;
        options.ti = 1;
        % compute coefficients
        MW = perform_wavelet_transform(M,Jmin, +1, options);

    case {'waveortho' 'wavortho'}
        amin = 2.5; amax = 3.6;
        options.ti = 0;
        % compute coefficients
        MW = perform_wavelet_transform(M, Jmin, +1, options);

    case 'curvelet'
        use_normalized = 0;
        niter = 15;
        if use_normalized
            amin = 2.2; % amin = 1.5;
            amax = 3.2;
        else
            amin = .5; amax = 1.2; % because of redundancy (no normalization of curvelets atoms)
        end

        if use_normalized
            F = ones(n);
            X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
            C = perform_curvelet_transform(X,options);

            % Compute norm of curvelets (exact)
            E = cell(size(C));
            for s=1:length(C)
                E{s} = cell(size(C{s}));
                for w=1:length(C{s})
                    A = C{s}{w};
                    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
                end
            end
        end

        % Take curvelet transform
        C = perform_curvelet_transform(M,options);

    case {'bandlet' 'bandlet-qt'}
        if not(isfield(options, 'G'))
            error('You should provide options.G');
        end
        G = options.G;
        niter = 15;
        amin = 2; amax = 3.5;
        amax = 4;
        if strcmp(type, 'bandlet-qt')
            bandcb = @perform_quadtree_transform;
        else
            bandcb = @perform_bandlet_transform;
        end
        % bandlets
        options.dir = 1;
        MB = feval(bandcb, M, G, +1, options );
    case {'gaussian'}
        % gaussian blurring
        amin = 1/sigma; amax = 10/sigma;
        thresh = '';
    otherwise 
        error('Unknown method');
end

if strcmp(thresh, 'soft')
    % lower the thresholds
    amin = amin/2; amax = amax/2;
end
if strcmp(thresh, 'block') || strcmp(thresh, 'block-soft')
    % lower the thresholds
    amin = amin/3; amax = amax/2;
end
if strcmp(thresh, 'block-hard')
    % lower the thresholds
    amin = amin*0; amax = amax*.4;
end

%% set up the list of threshold
Tlist = linspace(amin,amax,niter)*sigma;
Tlist = getoptions(options, 'Tlist', Tlist);
niter = length(Tlist);

%% perform the SNR optimization
M1 = {}; err = [];
for i = 1:niter
    T = Tlist(i);
    if verb
        progressbar(i,niter);
    end
    switch lower(type)
        case {'wavelet' 'waveti' 'wavti' 'waveortho' 'wavortho'}
            % wavelets
            MWT = perform_thresholding(MW, T, thresh, options);
            if not(thresh_low)
                % re-inject low frequencies
                if not(iscell(MWT))
                    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
                else
                    MWT{end} = MW{end};
                end
            end
            M1{end+1} = perform_wavelet_transform(MWT,Jmin,-1,options);
        case {'bandlet' 'bandlet-qt'}
            % perform thresholding
            MBT = perform_thresholding(MB, T, thresh);
            % reconstruct
            options.dir = -1;
            M1{end+1} = feval( bandcb, MBT, G, -1, options );
        case 'curvelet'
            % Apply thresholding
            if use_normalized
                Ct = C;
                for s = 2:length(C)
                    thr = T + sigma*(s == length(C));
                    for w = 1:length(C{s})
                        Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thr*E{s}{w});
                    end
                end
            else
                Ct = perform_thresholding(C,T,'hard');
            end
            % Take inverse curvelet transform
            M1{end+1} = perform_curvelet_transform(Ct,options);
            %            M1{end+1} = real( ifdct_wrapping(Ct, isreal, n,n) );
        case 'gaussian'          
            M1{end+1} = perform_blurring(M, T);
    end
    err(i) = snr(M0,M1{end});
end

% clf; plot(Tlist/sigma,err); axis tight;

% select optimal value
[pval,I] = max(err);
M1 = M1{I};
T = Tlist(I);

if I==1 || I==niter
    warning('Out of bound reached.');
end