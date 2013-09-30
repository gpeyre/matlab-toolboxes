function ML = perform_mca(M, components, options)

% perform_mca - perform MCA decomposition
%
%   ML = perform_mca(M, components, options);
%   
%   ML(:,:,i) is the layer optained by sparse decomposition in 
%       dictionary Di described by components{i}.
%
%   components is a cell array of structure.
%   cpt = components{i) discribe the ith transform used. 
%
%   It must contains: 
%       cpt.callback: a function of the type 
%           y = function(x,dir,options);
%       for instance you can set cpt.callback=@callback_atrou
%
%   It can contains:
%       cpt.options: the options given to the callback
%       cpt.threshold_factor: an amplification factor for the thresholding.
%       cpt.tv_correction and cpt.tv_weight: the MCA will perform an
%           additional TV minimization (usefull for the cartoon layer).
%
%   options contains options for the MCA process such as:
%       options.niter: number of iterations.
%       options.saverep: directory to save iterations
%       options.Tmax and options.Tmin: maximum and minimum threshold
%           values. Set Tmin=0 if there is no noise in the problem
%           (decomposition will be exact i.e. M=sum(ML,3)
%           Set Tmax to approximately the maximum value of the
%           coefficients of all the transforms of M.
%
%   The references for the Morphological Components Analysis are
%   
%       J. Bobin, Y. Moudden, J.-L. Starck and M. Elad,  
%       "Morphological Diversity and Source Separation",  
%       IEEE Transaction on Signal Processing , Vol 13, 7, pp 409--412, 2006.
%
%       M. Elad, J.-L Starck, D. Donoho and P. Querre, 
%       "Simultaneous Cartoon and Texture Image Inpainting using Morphological Component Analysis (MCA)", 
%       Journal on Applied and Computational Harmonic Analysis ACHA , Vol. 19, pp. 340-358, November 2005.
%
%   	J.-L. Starck, M. Elad, and D.L. Donoho, 
%       "Image Decomposition Via the Combination of Sparse Representation and a Variational Approach", 
%       IEEE Transaction on Image Processing , 14, 10, 2005.
%
%       J.-L. Starck, M. Elad, and D.L. Donoho, 
%       "Redundant Multiscale Transforms and their Application for Morphological Component Analysis", 
%       Advances in Imaging and Electron Physics , 132, 2004. 
%
%   See also: perform_iterative_thresholding.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 20;
end
if isfield(options, 'Tmax')
    Tmax = options.Tmax;
else
    Tmax = 0.2*max(abs(M(:)));
end
if isfield(options, 'Tmin')
    Tmin = options.Tmin;
else
    Tmin = 0;
end
if isfield(options, 'threshold')
    threshold = options.threshold;
else
    threshold = 'hard';
end
if isfield(options, 'saverep')
    saverep = options.saverep;
else
    saverep = [];
end
if isfield(options, 'InpaintingMask')
    InpaintingMask = options.InpaintingMask;
else
    InpaintingMask = [];
end

s = length(components);


% init the residual
n = size(M,1);
ML = zeros(n,n,s);
if not(isempty(InpaintingMask))
    q = sum(InpaintingMask(:)==Inf);
    for k=1:s
        Mk = ML(:,:,k);
        Mk(InpaintingMask==Inf) = rand( q,1 ) * max( M(InpaintingMask~=Inf) );
        ML(:,:,k) = Mk;
    end
end


MLsave = [];
for i=1:niter
    progressbar(i,niter);
    % current threshold
    T0 = Tmax + (i-1)/(niter-1)*(Tmin-Tmax);
    for k=1:s
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % current component
        cpt = components{k};
        % callback function
        callback = cpt.callback;
        % options
        clear opt; opt.null = 0;
        if isfield(cpt, 'options')
            opt = cpt.options;
        end
        % threshold factor amplification
        threshold_factor = 1;
        if isfield(cpt, 'threshold_factor')
            threshold_factor = cpt.threshold_factor;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute residual
        T = T0 * threshold_factor;
        opt.T = T;
        % compute residual
        sel = 1:s; sel(k) = [];
        R = M - sum( ML(:,:,sel), 3);
        if not(isempty(InpaintingMask))
            % enforce known values only outside
            Mk = ML(:,:,k);
            R(InpaintingMask==Inf) = Mk(InpaintingMask==Inf);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform forward transform
        RW = feval(callback, R,+1,opt);
        % perform thresholding
        if strcmp(cpt.name, 'wav')
            % do not threshold the low frequency
        	RW1 = perform_thresholding({RW{1:end-1}},T, threshold);
            RW = {RW1{:}, RW{end}};
        else
        	RW = perform_thresholding(RW,T, threshold);
        end
        % perform backward transform
        ML(:,:,k) = feval(callback, RW,-1,opt);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TV additional sparsening
        if isfield(cpt, 'tv_correction') && cpt.tv_correction==1
            if isfield(options, 'tv_weight')
                tv_weight = options.tv_weight;
            else
                tv_weight = 2.5*max(abs(M(:)))/256;
            end
            ML(:,:,k) = perform_tv_correction(ML(:,:,k),tv_weight);
            end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save result
        if not(isempty(saverep))
            MLsave(:,:,:,i) = ML;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save all the data
if not(isempty(saverep))
    % normalize the evolution
    vmax = max(max(max(MLsave,[],1),[],2),[],4); 
    vmax = repmat(vmax, [n n 1 niter]);
    vmin = min(min(min(MLsave,[],1),[],2),[],4); 
    vmin = repmat(vmin, [n n 1 niter]);
    MLsave = (MLsave-vmin) ./ (vmax-vmin);
    
    for k=1:s
        cpt = components{k};
        for i=1:niter
            if isfield(cpt, 'name')
                str = cpt.name;
            else
                str = ['layer' num2str(k)];
            end
            str = [str '-iter' num2string_fixeddigit(i,2)];
            warning off;
            imwrite( MLsave(:,:,k,i), [saverep str '.png'], 'png' );
            warning on;
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_hard_tresholding(x,t)

if iscell(x)
    for i=1:length(x)
        y{i} = perform_hard_tresholding(x{i},t);
    end
    return;
end

y = x .* (abs(x) > t);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_soft_tresholding(x,t)

if iscell(x)
    for i=1:length(x)
        y{i} = perform_soft_tresholding(x{i},t);
    end
    return;
end

s = abs(x) - t;
s = (s + abs(s))/2;
y = sign(x).*s;
