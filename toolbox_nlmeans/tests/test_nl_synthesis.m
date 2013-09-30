% test for texture synthesis using NL means
%
%   Copyright (c) 2007 Gabriel Peyre

path(path,'toolbox/');
path(path,'images/');
path(path,'../../images/');


if not(exist('name'))
% Color
% B&W
name = 'corral';
name = 'mures';

name = 'olives';
name = 'parrot';
name = 'pointilliste';
name = 'grass';
name = 'reptilskin';
name = 'fabric';
name = 'deadleaf';
end

%% test type
test = 'denoise';
test = 'inpainting';
test = 'synthesis';

if not(exist('T'))
%% smoothing factor
T = 0.15;
T = 0.06;
T = 1e-9; % strict recopy
end
% we will decay the size of the squares during synthesis
if not(exist('k_schedule'))
    k_schedule = [4];
    k_schedule = [4 3 2];
    k_schedule = [3 2];
    k_schedule = [5];
end

if not(exist('do_patchwise'))
    %% use averaging of overlapping patches
    do_patchwise = 1;
end

if T<1e-6
	Tstr = 'inf';
else
	Tstr = num2str(T);
end

repbase = ['results/' test '/' name '/'];
if not(exist(repbase))
    mkdir(repbase);
end
if length(k_schedule)==1
    repimg = [repbase 'k' num2str(k_schedule) '/'];
else   
    repimg = [repbase 'k' num2str(k_schedule(1)) '-' num2str(k_schedule(end)) '/'];
end
if not(exist(repimg))
    mkdir(repimg);
end
if strcmp(test, 'synthesis')
repquilt = [repbase 'quilting/'];
if not(exist(repquilt))
    mkdir(repquilt);
end
end

clear options;
% size of synthesized image
n = 256; m = n;
% size of original image
na = n;
na = 128;

n0 = [];
switch lower(name)
    case 'corral'
        n0 = 150;
    case 'dunes'
        n0 = 200;
    case 'warped_grid'
        n0 = 450;
    case 'wood'
        n0 = 230;
    case 'grass-big'
        n0 = 200;
    case 'dead_leaf'
        n0 = 400;
    case 'fabric'
        n0 = 400;
    case 'simoncelli3'
        n0 = 300;
    case 'simoncelli7'
        n0 = 180;
    case 'text'
        n0 = 200;
    case 'turbulence1'
        n0 = 150;
    otherwise
        n0 = [];
end

%% load the image
if not(strcmp(name, 'parrot'))
Ma = load_image(name, n0);
na = min(size(Ma,1),na); ma = na;   % size of the original image
Ma = rescale( crop(Ma, [ma na]) );
s = size(Ma,3); % number of colors
if not(strcmp(test, 'synthesis'))
    n = na;
end
end

%% parameters for nl-means
% we use a very low variance to perform strict recopy
options.T = T; 
if strcmp(test, 'synthesis')
    options.max_dist = 12;  % distance for search
else
    options.max_dist = 15;  % distance for search
end
options.ndims = 15; % number of dimension for matching (PCA) 30 is fine


% perform or not wavelet histogram matching
use_wavelet_histo = 1;
use_spacial_histo = 1;  
if s>1
    use_spacial_histo = 0;
end

disp(['--> Synthesizing ' name ', ave=' num2str(do_patchwise) ', T=' Tstr '.']);

if strcmp(test, 'inpainting')
    % determine the mask
    if strcmp(name, 'deadleaf') || strcmp(name, 'fabric') || strcmp(name, 'grass') || strcmp(name, 'reptilskin')
        path(path, 'images/inpainting-perez/');
        Ma = load_image( [name '-masked']); Ma = rescale(sum(Ma,3));
        na = size(Ma,1); n = na;
        mask_copy = double(Ma==0);
    elseif strcmp(name, 'parrot')
        path(path, 'images/inpainting-fadili/');
        Ma = load_image( name); 
        c = round(size(Ma)*0.38); c(1)=c(1)-35;
        Ma = rescale( crop(Ma,na, c ) );
        n = na; s = size(Ma,3);
        mask_copy = crop( load_image( [name '-mask']), na, c); 
        mask_copy = rescale(-sum(mask_copy,3));
        mask_copy = double(mask_copy==1);        
    else
        if not(exist('mask_copy'))
            options.r = 5; options.mode = 'line';
            mask_copy = grab_inpainting_mask(Ma, options); mask_copy=double(sum(mask_copy,3)==Inf);
        end
    end
    use_wavelet_histo = 0;
    use_spacial_histo = 0;
    I0 = find(mask_copy<0.5); I1 = find(mask_copy>=0.5);
    Irecopy = []; Iremove = [];
    for j=1:s
        Irecopy = [Irecopy; I0+(j-1)*n^2];
        Iremove = [Iremove; I1+(j-1)*n^2];
    end
    % remove pixels
    M0 = Ma;
    Ma(Iremove) = 1;
    warning off;
    imwrite(clamp(M0), [repbase name '-target.png'], 'png');
    imwrite(clamp(Ma), [repbase name '-mask.png'], 'png');
    warning on;
    Ma(Iremove) = rand(length(Iremove),1);
    % compute mask for inpainting
    mask_process = double( perform_convolution(mask_copy,ones(3))>0 );
    options.mask_process = mask_process;
    options.mask_copy = mask_copy;
    options.exclude_self = 1;
end

options.wmax = max(k_schedule);
% initialization for windows search, just random
if strcmp(test, 'synthesis')
	M = randn(m,n,s);
	options.Vx = floor( rand(m,n)*(ma-1) ) + 1;
	options.Vy = floor( rand(m,n)*(na-1) ) + 1;
else
	M = Ma;
end

% do a little diffusion as initialization
niter_diffusion = 50; sigma = 1;
if 0
if strcmp(test, 'inpainting')
    for i=1:niter_diffusion
        M = perform_blurring(M,sigma);
        M(Irecopy) = Ma(Irecopy);
    end
end
end
    
    
%% save originals
warning off;
imwrite(clamp(Ma), [repbase name '-original.png'], 'png');
imwrite(clamp(M), [repimg name '-synth-00.png'], 'png');
warning on;

%% initialization
if use_wavelet_histo
    options.dotransform = 1; Jmin = 4;
    M = perform_wavelet_matching(M,Ma,options);
end
if use_spacial_histo
    M = perform_histogram_equalization(M, Ma, options);
end

% number of iterations for synthesis
niter = 12;
options.do_patchwise = do_patchwise; % use or not local averaging during recopy

%% iterations
for i=1:niter
    progressbar(i,niter);
    % set windows half width
    ik = floor( (i-1)/niter*length(k_schedule) )+1;
    options.k = k_schedule(ik);
    options.Ma = Ma;
    if strcmp(test, 'inpainting')
        options.Ma = M;
    end
    [M,Vx,Vy,new_mask] = perform_nl_means(M, options);
    if i>100 && mod(i,6)==0
        % update the mask
        options.mask_copy = 1 - double( perform_convolution(1-options.mask_copy,ones(3))>0 );
    end
    if strcmp(test, 'synthesis')
		options.Vx = Vx; options.Vy = Vy;
    end
    if strcmp(test, 'inpainting')
        M(Irecopy) = Ma(Irecopy);
    end
    if use_wavelet_histo
        M = perform_wavelet_matching(M,Ma,options);
    end 
    if use_spacial_histo
        M = perform_histogram_equalization(M, Ma, options);
    end
    M = clamp(M);
    % save result
    warning off;
    imwrite(clamp(M), [repimg name '-synth-' num2string_fixeddigit(i,2) '.png'], 'png');
    warning on;
end


warning off;
imwrite(clamp(M), [repimg name '-synthesis-nlmeans.png'], 'png');
warning on;


%% display final result
n1 = na;
clf;
subplot(1,2,1);
d = (n1-n)/2;
imagesc(Ma); axis image; axis([1-d n1-d 1-d n1-d]); axis off;
title('Original');
subplot(1,2,2);
imagesc(M); axis image; axis off;
title('Synthetised');
if size(M,3)==1
    colormap gray(256);
end

saveas(gcf, [repbase name '-aver' num2str(do_patchwise) '-' test '-' Tstr '.png'], 'png');

if not(exist('do_quilting'))
    do_quilting = 1;
end

if strcmp(test, 'synthesis') && do_quilting
    % perform a traditional texture synthesis with quilting
    for tilesize = 12 % [8 16 20]
        overlap = max(3,ceil(0.25*tilesize));
        ntiles = ceil((n-overlap)/(tilesize-overlap));
        disp(['--> Synthesis using quilting, w=' num2str(tilesize) '.']);
        clf;
        Mquilt = perform_synthesis_quilting(Ma, tilesize, ntiles, overlap);
        Mquilt = crop(Mquilt, n);
        warning off
        imwrite(clamp(Mquilt), [repquilt name '-synthesis-quilting-w' num2str(tilesize) '.png'], 'png');
        warning on;
    end;
end