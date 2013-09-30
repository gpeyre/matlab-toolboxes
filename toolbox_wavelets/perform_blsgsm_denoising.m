function y = perform_blsgsm_denoising(x, options)

% perform_blsgsm_denoising - denoise an image using BLS-GSM
%
% y = perform_blsgsm_denoising(x, options);
%
%   BLS-GSM stands for "Bayesian Least Squares - Gaussian Scale Mixture". 
%
%   This function is a wrapper for the code of J.Portilla.
%
%   You can change the following optional parameters
%   options.Nor: Number of orientations (for X-Y separable wavelets it can
%       only be 3)
%   options.repres1: Type of pyramid ('uw': shift-invariant version of an orthogonal wavelet, in this case)
%   options.repres2: Type of wavelet (daubechies wavelet, order 2N, for 'daubN'; in this case, 'Haar')
%
%   Other parameter that should not be changed unless really needed.
%   options.blSize	    	n x n coefficient neighborhood (n must be odd): 
%   options.parent			including or not (1/0) in the
%       neighborhood a coefficient from the same spatial location
%   options.boundary        Boundary mirror extension, to avoid boundary artifacts 
%   options.covariance     	Full covariance matrix (1) or only diagonal elements (0).
%   options.optim          	Bayes Least Squares solution (1), or MAP-Wiener solution in two steps (0)
%
%   Default parameters favors speed.
%   To get the optimal results, one should use:
%       options.Nor = 8;                8 orientations
%       options.repres1 = 'fs';         Full Steerable Pyramid, 5 scales for 512x512
%       options.repres2 = '';           Dummy parameter when using repres1 = 'fs'   
%       options.parent = 1;             Include a parent in the neighborhood
%
%   J Portilla, V Strela, M Wainwright, E P Simoncelli, 
%   "Image Denoising using Scale Mixtures of Gaussians in the Wavelet Domain,"
%   IEEE Transactions on Image Processing. vol 12, no. 11, pp. 1338-1351,
%   November 2003 
%
%   Javier Portilla, Universidad de Granada, Spain
%   Adapted by Gabriel Peyre in 2007

options.null = 0;

if isfield(options, 'sigma')
    sig = options.sigma;
else
    sig = perform_noise_estimation(x)
end


[Ny,Nx] = size(x);
% Pyramidal representation parameters
Nsc = ceil(log2(min(Ny,Nx)) - 4);  % Number of scales (adapted to the image size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
if isfield(options, 'Nor')
    Nor = options.Nor;
else
    Nor = 3;				            % Number of orientations (for X-Y separable wavelets it can only be 3)
end
if isfield(options, 'repres1')
    repres1 = options.repres1;
else
    repres1 = 'uw';                     % Type of pyramid (shift-invariant version of an orthogonal wavelet, in this case)
end

if isfield(options, 'repres2')
    repres2 = options.repres2;
else
    repres2 = 'daub1';                  % Type of wavelet (daubechies wavelet, order 2N, for 'daubN'; in this case, 'Haar')
end

% Model parameters (optimized: do not change them unless you are an advanced user with a deep understanding of the theory)
if isfield(options, 'blSize')
    blSize = options.blSize;
else
    blSize = [3 3];	    % n x n coefficient neighborhood of spatial neighbors within the same subband
end

if isfield(options, 'parent')
    parent = options.parent;
else
    parent = 0;			% including or not (1/0) in the neighborhood a coefficient from the same spatial location
end
if isfield(options, 'boundary')
    boundary = options.boundary;
else
    boundary = 1;		% Boundary mirror extension, to avoid boundary artifacts
end
if isfield(options, 'covariance')
    covariance = options.covariance;
else
    covariance = 1;     % Full covariance matrix (1) or only diagonal elements (0).
end
if isfield(options, 'optim')
    optim = options.optim;
else
    optim = 1;          % Bayes Least Squares solution (1), or MAP-Wiener solution in two steps (0)
end

PS = ones(size(x));	% power spectral density (in this case, flat, i.e., white noise)
seed = 0;

% Uncomment the following 4 code lines for reproducing the results of our IEEE Trans. on Im. Proc., Nov. 2003 paper
% This configuration is slower than the previous one, but it gives slightly better results (SNR)
% on average for the test images "lena", "barbara", and "boats" used in the
% cited article.
% Nor = 8;                           % 8 orientations
% repres1 = 'fs';                    % Full Steerable Pyramid, 5 scales for 512x512
% repres2 = '';                      % Dummy parameter when using repres1 = 'fs'   
% parent = 1;                        % Include a parent in the neighborhood


% Call the denoising function
y = denoi_BLS_GSM(x, sig, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed);




function im_d = denoi_BLS_GSM(im, sig, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed);

% [im_d,im,SNR_N,SNR,PSNR] = denoi_BLS_GSM(im, sig, ft, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed);
%
%	im:	input noisy image
%	sig:	standard deviation of noise
%	PS:	Power Spectral Density of noise ( fft2(autocorrelation) )
%			NOTE: scale factors do not matter. Default is white.
%	blSize: 2x1 or 1x2 vector indicating local neighborhood
%		([sY sX], default is [3 3])
%   parent: Use parent yes/no (1/0)
%	Nsc:	Number of scales
%   Nor:  Number of orientations. For separable wavelets this MUST be 3.
%   covariance: Include / Not Include covariance in the GSM model (1/0)
%   optim: BLS / MAP-Wiener(2-step) (1/0)
%   repres1: Possible choices for representation:
%           'w':    orthogonal wavelet
%                   (uses buildWpyr, reconWpyr)
%                   repres2 (optional):
%                      haar:              - Haar wavelet.
%                      qmf8, qmf12, qmf16 - Symmetric Quadrature Mirror Filters [Johnston80]
%                      daub2,daub3,daub4  - Daubechies wavelet [Daubechies88] (#coef = 2N, para daubN).
%                      qmf5, qmf9, qmf13: - Symmetric Quadrature Mirror Filters [Simoncelli88,Simoncelli90]
%           'uw':   undecimated orthogonal wavelet, Daubechies, pyramidal version
%                   (uses buildWUpyr, reconWUpyr).
%                   repres2 (optional): 'daub<N>', where <N> is a positive integer (e.g., 2)
%           's':    steerable pyramid [Simoncelli&Freeman95].
%                   (uses buildSFpyr, reconSFpyr)
%           'fs':   full steerable pyramid [Portilla&Simoncelli02].
%                   (uses buildFullSFpyr2, reconsFullSFpyr2)
%   seed (optional):    Seed used for generating the Gaussian noise (when ft == 0)
%                       By default is 0.
%
%   im_d: denoising result

% Javier Portilla, Univ. de Granada, 5/02
% revision 31/03/2003
% revision 7/01/2004
% Last revision 15/11/2004

if ~exist('blSize'),
    blSzX = 3;	% Block size
    blSzY = 3;
else
    blSzY = blSize(1);
    blSzX = blSize(2);
end

if (blSzY/2==floor(blSzY/2))|(blSzX/2==floor(blSzX/2)),
    error('Spatial dimensions of neighborhood must be odd!');
end

if ~exist('PS'),
    no_white = 0;   % Power spectral density of noise. Default is white noise
else
    no_white = 1;
end

if ~exist('parent'),
    parent = 1;
end

if ~exist('boundary'),
    boundary = 1;
end

if ~exist('Nsc'),
    Nsc = 4;
end

if ~exist('Nor'),
    Nor = 8;
end

if ~exist('covariance'),
    covariance = 1;
end

if ~exist('optim'),
    optim = 1;
end

if ~exist('repres1'),
    repres1 = 'fs';
end

if ~exist('repres2'),
    repres2 = '';
end

if (((repres1=='w') | (repres1=='uw')) & (Nor~=3)),
    warning('For X-Y separable representations Nor must be 3. Nor = 3 assumed.');
    Nor = 3;
end

if ~exist('seed'),
    seed = 0;
end

[Ny Nx] = size(im);

% We ensure that the processed image has dimensions that are integer
% multiples of 2^(Nsc+1), so it will not crash when applying the
% pyramidal representation. The idea is padding with mirror reflected
% pixels (thanks to Jesus Malo for this idea).

Npy = ceil(Ny/2^(Nsc+1))*2^(Nsc+1);
Npx = ceil(Nx/2^(Nsc+1))*2^(Nsc+1);

if Npy~=Ny | Npx~=Nx,
    Bpy = Npy-Ny;
    Bpx = Npx-Nx;
    im = bound_extension(im,Bpy,Bpx,'mirror');
    im = im(Bpy+1:end,Bpx+1:end);	% add stripes only up and right
end


% size of the extension for boundary handling
if (repres1 == 's') | (repres1 == 'fs'),
    By = (blSzY-1)*2^(Nsc-2);
    Bx = (blSzX-1)*2^(Nsc-2);
else
    By = (blSzY-1)*2^(Nsc-1);
    Bx = (blSzX-1)*2^(Nsc-1);
end

if ~no_white,       % White noise
    PS = ones(size(im));
end

% As the dimensions of the power spectral density (PS) support and that of the
% image (im) do not need to be the same, we have to adapt the first to the
% second (zero padding and/or cropping).

PS = fftshift(PS);
isoddPS_y = (size(PS,1)~=2*(floor(size(PS,1)/2)));
isoddPS_x = (size(PS,2)~=2*(floor(size(PS,2)/2)));
PS = PS(1:end-isoddPS_y, 1:end-isoddPS_x);          % ensures even dimensions for the power spectrum
PS = fftshift(PS);

[Ndy,Ndx] = size(PS);   % dimensions are even

delta = real(ifft2(sqrt(PS)));
delta = fftshift(delta);
aux = delta;
delta = zeros(Npy,Npx);
if (Ndy<=Npy)&(Ndx<=Npx),
    delta(Npy/2+1-Ndy/2:Npy/2+Ndy/2,Npx/2+1-Ndx/2:Npx/2+Ndx/2) = aux;
elseif (Ndy>Npy)&(Ndx>Npx),
    delta = aux(Ndy/2+1-Npy/2:Ndy/2+Npy/2,Ndx/2+1-Npx/2:Ndx/2+Npx/2);
elseif (Ndy<=Npy)&(Ndx>Npx),
    delta(Npy/2+1-Ndy/2:Npy/2+1+Ndy/2-1,:) = aux(:,Ndx/2+1-Npx/2:Ndx/2+Npx/2);
elseif (Ndy>Npy)&(Ndx<=Npx),
    delta(:,Npx/2+1-Ndx/2:Npx/2+1+Ndx/2-1) = aux(Ndy/2+1-Npy/2:Ndy/2+1+Npy/2-1,:);
end

if repres1 == 'w',
    PS = abs(fft2(delta)).^2;
    PS = fftshift(PS);
    % noise, to be used only with translation variant transforms (such as orthogonal wavelet)
    delta = real(ifft2(sqrt(PS).*exp(j*angle(fft2(randn(size(PS)))))));
end

%Boundary handling: it extends im and delta
if boundary,
    im = bound_extension(im,By,Bx,'mirror');
    if repres1 == 'w',
        delta = bound_extension(delta,By,Bx,'mirror');
    else
        aux = delta;
        delta = zeros(Npy + 2*By, Npx + 2*Bx);
        delta(By+1:By+Npy,Bx+1:Bx+Npx) = aux;
    end
else
    By=0;Bx=0;
end

delta = delta/sqrt(mean2(delta.^2));    % Normalize the energy (the noise variance is given by "sig")
delta = sig*delta;                      % Impose the desired variance to the noise


% main
t1 = clock;
if repres1 == 's',  % standard steerable pyramid
    im_d = decomp_reconst(im, Nsc, Nor, [blSzX blSzY], delta, parent,covariance,optim,sig);
elseif repres1 == 'fs', % full steerable pyramid
    im_d = decomp_reconst_full(im, Nsc, Nor, [blSzX blSzY], delta, parent, covariance, optim, sig);
elseif repres1 == 'w',  % orthogonal wavelet
    if ~exist('repres2'),
        repres2 = 'daub1';
    end
    filter = repres2;
    im_d = decomp_reconst_W(im, Nsc, filter, [blSzX blSzY], delta, parent, covariance, optim, sig);
elseif repres1 == 'uw',    % undecimated daubechies wavelet
    if ~exist('repres2'),
        repres2 = 'daub1';
    end
    if repres2(1:4) == 'haar',
        daub_order = 2;
    else
        daub_order = 2*str2num(repres2(5));
    end
    im_d = decomp_reconst_WU(im, Nsc, daub_order, [blSzX blSzY], delta, parent, covariance, optim, sig);
else
    error('Invalid representation parameter. See help info.');
end
t2 = clock;
elaps = t2 - t1;
elaps(4)*3600+elaps(5)*60+elaps(6); % elapsed time, in seconds

im_d = im_d(By+1:By+Npy,Bx+1:Bx+Npx);
im_d = im_d(1:Ny,1:Nx);

function fh = decomp_reconst_full(im,Nsc,Nor,block,noise,parent,covariance,optim,sig);

% Decompose image into subbands, denoise, and recompose again.
%		fh = decomp_reconst(im,Nsc,Nor,block,noise,parent,covariance,optim,sig);
%       covariance:	 are we considering covariance or just variance?
%       optim:		 for choosing between BLS-GSM (optim = 1) and MAP-GSM (optim = 0)
%       sig:        standard deviation (scalar for uniform noise or matrix for spatially varying noise)
% Version using the Full steerable pyramid (2) (High pass residual
% splitted into orientations).

% JPM, Univ. de Granada, 5/02
% Last Revision: 11/04

if (block(1)/2==floor(block(1)/2))|(block(2)/2==floor(block(2)/2)),
    error('Spatial dimensions of neighborhood must be odd!');
end

if ~exist('parent'),
    parent = 1;
end

if ~exist('covariance'),
    covariance = 1;
end

if ~exist('optim'),
    optim = 1;
end

if ~exist('sig'),
    sig = sqrt(mean(noise.^2));
end

[pyr,pind] = buildFullSFpyr2(im,Nsc,Nor-1);
[pyrN,pind] = buildFullSFpyr2(noise,Nsc,Nor-1);
pyrh = real(pyr);
Nband = size(pind,1)-1;
for nband = 2:Nband, % everything except the low-pass residual
    fprintf('%d % ',round(100*(nband-1)/(Nband-1)))
    aux = pyrBand(pyr, pind, nband);
    auxn = pyrBand(pyrN, pind, nband);
    [Nsy,Nsx] = size(aux);
    prnt = parent & (nband < Nband-Nor);   % has the subband a parent?
    BL = zeros(size(aux,1),size(aux,2),1 + prnt);
    BLn = zeros(size(aux,1),size(aux,2),1 + prnt);
    BL(:,:,1) = aux;
    BLn(:,:,1) = auxn*sqrt(((Nsy-2)*(Nsx-2))/(Nsy*Nsx));     % because we are discarding 2 coefficients on every dimension
    if prnt,
        aux = pyrBand(pyr, pind, nband+Nor);
        auxn = pyrBand(pyrN, pind, nband+Nor);
        if nband>Nor+1,     % resample 2x2 the parent if not in the high-pass oriented subbands.
            aux = real(expand(aux,2));
            auxn = real(expand(auxn,2));
        end
        BL(:,:,2) = aux;
        BLn(:,:,2) = auxn*sqrt(((Nsy-2)*(Nsx-2))/(Nsy*Nsx)); % because we are discarding 2 coefficients on every dimension
    end

    sy2 = mean2(BL(:,:,1).^2);
    sn2 = mean2(BLn(:,:,1).^2);
    if sy2>sn2,
        SNRin = 10*log10((sy2-sn2)/sn2);
    else
        disp('Signal is not detectable in noisy subband');
    end

    % main
    BL = denoi_BLS_GSM_band(BL,block,BLn,prnt,covariance,optim,sig);
    pyrh(pyrBandIndices(pind,nband)) = BL(:)';
end
fh = reconFullSFpyr2(pyrh,pind);


function fh = decomp_reconst_W(im,Nsc,filter,block,noise,parent,covariance,optim,sig);

% Decompose image into subbands, denoise using BLS-GSM method, and recompose again.
%		fh = decomp_reconst(im,Nsc,filter,block,noise,parent,covariance,optim,sig);
%       im:         image
%       Nsc:        number of scales
%       filter:     type of filter used (see namedFilters)
%       block:      2x1 vector indicating the dimensions (rows and columns) of the spatial neighborhood 
%       noise:      signal with the same autocorrelation as the noise
%       parent:     include (1) or not (0) a coefficient from the immediately coarser scale in the neighborhood
%       covariance:	 are we considering covariance or just variance?
%       optim:		 for choosing between BLS-GSM (optim = 1) and MAP-GSM (optim = 0)
%       sig:        standard deviation (scalar for uniform noise or matrix for spatially varying noise)
% Version using a critically sampled pyramid (orthogonal wavelet), as implemented in MatlabPyrTools (Eero).

% JPM, Univ. de Granada, 3/03

if (block(1)/2==floor(block(1)/2))|(block(2)/2==floor(block(2)/2)),
   error('Spatial dimensions of neighborhood must be odd!');
end   

if ~exist('parent'),
        parent = 1;
end

if ~exist('covariance'),
        covariance = 1;
end

if ~exist('optim'),
        optim = 1;
end

if ~exist('sig'),
        sig = sqrt(mean(noise.^2));
end

Nor = 3;    % number of orientations: vertical, horizontal and mixed diagonals (for compatibility)

[pyr,pind] = buildWpyr(im,Nsc,filter,'circular');
[pyrN,pind] = buildWpyr(noise,Nsc,filter,'circular');
pyrh = pyr;
Nband = size(pind,1);
for nband = 1:Nband-1, % everything except the low-pass residual
  fprintf('%d % ',round(100*(nband-1)/(Nband-1)))
  aux = pyrBand(pyr, pind, nband);
  auxn = pyrBand(pyrN, pind, nband);
  prnt = parent & (nband < Nband-Nor);   % has the subband a parent?
  BL = zeros(size(aux,1),size(aux,2),1 + prnt);
  BLn = zeros(size(aux,1),size(aux,2),1 + prnt);
  BL(:,:,1) = aux;
  BLn(:,:,1) = auxn;
  if prnt,
  	aux = pyrBand(pyr, pind, nband+Nor);
    auxn = pyrBand(pyrN, pind, nband+Nor);
    aux = real(expand(aux,2));
    auxn = real(expand(auxn,2));
    BL(:,:,2) = aux;
  	BLn(:,:,2) = auxn;
  end
  
  sy2 = mean2(BL(:,:,1).^2);
  sn2 = mean2(BLn(:,:,1).^2);
  if sy2>sn2,
     SNRin = 10*log10((sy2-sn2)/sn2);
  else
     disp('Signal is not detectable in noisy subband');
  end   
  
  % main
  BL = denoi_BLS_GSM_band(BL,block,BLn,prnt,covariance,optim,sig);
  pyrh(pyrBandIndices(pind,nband)) = BL(:)';
end
fh = reconWpyr(pyrh,pind,filter,'circular');


function fh = decomp_reconst_WU(im,Nsc,daub_order,block,noise,parent,covariance,optim,sig);

% Decompose image into subbands (undecimated wavelet), denoise, and recompose again.
%		fh = decomp_reconst_wavelet(im,Nsc,daub_order,block,noise,parent,covariance,optim,sig);
%       im :         image
%       Nsc:         Number of scales
%       daub_order:  Order of the daubechie fucntion used (must be even).
%       block:       size of neighborhood within each undecimated subband.
%       noise:       image having the same autocorrelation as the noise (e.g., a delta, for white noise)
%       parent:      are we including the coefficient at the central location at the next coarser scale?
%       covariance:	 are we considering covariance or just variance?
%       optim:		 for choosing between BLS-GSM (optim = 1) and MAP-GSM (optim = 0)
%       sig:        standard deviation (scalar for uniform noise or matrix for spatially varying noise)

% Javier Portilla, Univ. de Granada, 3/03
% Revised: 11/04 

if (block(1)/2==floor(block(1)/2))|(block(2)/2==floor(block(2)/2)),
   error('Spatial dimensions of neighborhood must be odd!');
end   

if ~exist('parent'),
        parent = 1;
end

if ~exist('covariance'),
        covariance = 1;
end

if ~exist('optim'),
        optim = 1;
end

if ~exist('sig'),
        sig = sqrt(mean(noise.^2));
end

Nor = 3;    % Number of orientations: vertical, horizontal and (mixed) diagonals.

[pyr,pind] = buildWUpyr(im,Nsc,daub_order);
[pyrN,pind] = buildWUpyr(noise,Nsc,daub_order);
pyrh = real(pyr);
Nband = size(pind,1)-1;

for nband = 2:Nband, % everything except the low-pass residual
  % fprintf('%d % ',round(100*(nband-1)/(Nband-1)))
  progressbar(nband-1,Nband-1);
  aux = pyrBand(pyr, pind, nband);
  auxn = pyrBand(pyrN, pind, nband);
  [Nsy, Nsx] = size(aux);
  prnt = parent & (nband < Nband-Nor);   % has the subband a parent?
  BL = zeros(size(aux,1),size(aux,2),1 + prnt);
  BLn = zeros(size(aux,1),size(aux,2),1 + prnt);
  BL(:,:,1) = aux;
  BLn(:,:,1) = auxn*sqrt(((Nsy-2)*(Nsx-2))/(Nsy*Nsx));     % because we are discarding 2 coefficients on every dimension
  if prnt,
  	aux = pyrBand(pyr, pind, nband+Nor);
    auxn = pyrBand(pyrN, pind, nband+Nor);
    if nband>Nor+1,     % resample 2x2 the parent if not in the high-pass oriented subbands.
	   aux = real(expand(aux,2));
       auxn = real(expand(auxn,2));
    end    
  	BL(:,:,2) = aux;
    BLn(:,:,2) = auxn*sqrt(((Nsy-2)*(Nsx-2))/(Nsy*Nsx)); % because we are discarding 2 coefficients on every dimension   
  end
  
  sy2 = mean2(BL(:,:,1).^2);
  sn2 = mean2(BLn(:,:,1).^2);
  if sy2>sn2,
     SNRin = 10*log10((sy2-sn2)/sn2);
  else
     disp('Signal is not detectable in noisy subband');
  end   
  
  % main
  BL = denoi_BLS_GSM_band(BL,block,BLn,prnt,covariance,optim,sig);
  pyrh(pyrBandIndices(pind,nband)) = BL(:)';
end
fh = reconWUpyr(pyrh,pind,daub_order);


function fh = decomp_reconst(fn,Nsc,Nor,block,noise,parent,covariance,optim,sig);

% Decompose image into subbands, denoise, and recompose again.
%	fh = decomp_reconst(fn,Nsc,Nor,block,noise,parent);

% Javier Portilla, Univ. de Granada, 5/02
% Last revision: 11/04

if (block(1)/2==floor(block(1)/2))|(block(2)/2==floor(block(2)/2)),
   error('Spatial dimensions of neighborhood must be odd!');
end   

if ~exist('parent'),
        parent = 1;
end

if ~exist('covariance'),
        covariance = 1;
end

if ~exist('optim'),
        optim = 1;
end

[pyr,pind] = buildSFpyr(fn,Nsc,Nor-1);
[pyrN,pind] = buildSFpyr(noise,Nsc,Nor-1);
pyrh = real(pyr);
Nband = size(pind,1);
for nband = 1:Nband -1,
  fprintf('%d % ',round(100*nband/(Nband-1)))
  aux = pyrBand(pyr, pind, nband);
  auxn = pyrBand(pyrN, pind, nband);
  [Nsy,Nsx] = size(aux);  
  prnt = parent & (nband < Nband-1-Nor) & (nband>1);
  BL = zeros(size(aux,1),size(aux,2),1 + prnt);
  BLn = zeros(size(aux,1),size(aux,2),1 + prnt);
  BL(:,:,1) = aux;
  BLn(:,:,1) = auxn*sqrt(((Nsy-2)*(Nsx-2))/(Nsy*Nsx));     % because we are discarding 2 coefficients on every dimension  
  if prnt,
  	aux = pyrBand(pyr, pind, nband+Nor);
	aux = real(expand(aux,2));
  	auxn = pyrBand(pyrN, pind, nband+Nor);
	auxn = real(expand(auxn,2));
  	BL(:,:,2) = aux;
    BLn(:,:,2) = auxn*sqrt(((Nsy-2)*(Nsx-2))/(Nsy*Nsx)); % because we are discarding 2 coefficients on every dimension       
  end
  
  sy2 = mean2(BL(:,:,1).^2);
  sn2 = mean2(BLn(:,:,1).^2);
  if sy2>sn2,
     SNRin = 10*log10((sy2-sn2)/sn2);
  else
     disp('Signal is not detectable in noisy subband');
  end
  
  % main
  BL = denoi_BLS_GSM_band(BL,block,BLn,prnt,covariance,optim,sig);
  pyrh(pyrBandIndices(pind,nband)) = BL(:)';
end
fh = reconSFpyr(pyrh,pind);


function x_hat = denoi_BLS_GSM_band(y,block,noise,prnt,covariance,optim,sig);

% It solves for the BLS global optimum solution, using a flat (pseudo)prior for log(z)
% 		  x_hat = denoi_BLS_GSM_band(y,block,noise,prnt,covariance,optim,sig);
%
%       prnt:  Include/ Not Include parent (1/0)
%       covariance: Include / Not Include covariance in the GSM model (1/0)
%       optim: BLS / MAP-Wiener(2-step) (1/0)

% JPM, Univ. de Granada, 5/02
% Last revision: JPM, 4/03


if ~exist('covariance'),
        covariance = 1;
end

if ~exist('optim'),
        optim = 1;
end

[nv,nh,nb] = size(y);

nblv = nv-block(1)+1;	% Discard the outer coefficients 
nblh = nh-block(2)+1;   % for the reference (centrral) coefficients (to avoid boundary effects)
nexp = nblv*nblh;			% number of coefficients considered
zM = zeros(nv,nh);		% hidden variable z
x_hat = zeros(nv,nh);	% coefficient estimation
N = prod(block) + prnt; % size of the neighborhood

Ly = (block(1)-1)/2;		% block(1) and block(2) must be odd!
Lx = (block(2)-1)/2;
if (Ly~=floor(Ly))|(Lx~=floor(Lx)),
   error('Spatial dimensions of neighborhood must be odd!');
end   
cent = floor((prod(block)+1)/2);	% reference coefficient in the neighborhood 
                                 % (central coef in the fine band)

Y = zeros(nexp,N);		% It will be the observed signal (rearranged in nexp neighborhoods)
W = zeros(nexp,N);		% It will be a signal with the same autocorrelation as the noise

foo = zeros(nexp,N);

% Compute covariance of noise from 'noise'
n = 0;
for ny = -Ly:Ly,	% spatial neighbors
	for nx = -Lx:Lx,
		n = n + 1;
		foo = shift(noise(:,:,1),[ny nx]);
		foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
		W(:,n) = vector(foo);
	end
end
if prnt,	% parent
	n = n + 1;
	foo = noise(:,:,2);
	foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
	W(:,n) = vector(foo);
end

C_w = innerProd(W)/nexp;
sig2 = mean(diag(C_w(1:N-prnt,1:N-prnt)));	% noise variance in the (fine) subband

clear W;
if ~covariance,
   if prnt,
        C_w = diag([sig2*ones(N-prnt,1);C_w(N,N)]);
   else
        C_w = diag(sig2*ones(N,1));
   end
end    


% Rearrange observed samples in 'nexp' neighborhoods 
n = 0;
for ny=-Ly:Ly,	% spatial neighbors
	for nx=-Lx:Lx,
		n = n + 1;
		foo = shift(y(:,:,1),[ny nx]);
		foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
		Y(:,n) = vector(foo);
	end
end
if prnt,	% parent
	n = n + 1;
	foo = y(:,:,2);
	foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
	Y(:,n) = vector(foo);
end
clear foo

% For modulating the local stdv of noise
if exist('sig') & prod(size(sig))>1,
    sig = max(sig,sqrt(1/12));   % Minimum stdv in quantified (integer) pixels
    subSampleFactor = log2(sqrt(prod(size(sig))/(nv*nh)));
    zW = blurDn(reshape(sig, size(zM)*2^subSampleFactor)/2^subSampleFactor,subSampleFactor);
    zW = zW.^2;
    zW = zW/mean2(zW); % Expectation{zW} = 1
    z_w = vector(zW(Ly+1:Ly+nblv,Lx+1:Lx+nblh));
end    

[S,dd] = eig(C_w);
S = S*real(sqrt(dd));	% S*S' = C_w
iS = pinv(S);
clear noise

C_y = innerProd(Y)/nexp;
sy2 = mean(diag(C_y(1:N-prnt,1:N-prnt))); % observed (signal + noise) variance in the subband
C_x = C_y - C_w;			% as signal and noise are assumed to be independent
[Q,L] = eig(C_x);
% correct possible negative eigenvalues, without changing the overall variance
L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
C_x = Q*L*Q';
   
sx2 = sy2 - sig2;			% estimated signal variance in the subband
sx2 = sx2.*(sx2>0); % + (sx2<=0); 
if ~covariance,
   if prnt,
        C_x = diag([sx2*ones(N-prnt,1);C_x(N,N)]);
   else
        C_x = diag(sx2*ones(N,1));
   end
end    
[Q,L] = eig(iS*C_x*iS');	 	% Double diagonalization of signal and noise
la = diag(L);						% eigenvalues: energy in the new represetnation.
la = real(la).*(real(la)>0);

% Linearly transform the observations, and keep the quadratic values (we do not model phase)

V = Q'*iS*Y';
clear Y;
V2 = (V.^2).';
M = S*Q;
m = M(cent,:);


% Compute p(Y|log(z))

if 1,   % non-informative prior
    lzmin = -20.5;
    lzmax = 3.5;
    step = 2;
else    % gamma prior for 1/z
    lzmin = -6;
    lzmax = 4;
    step = 0.5;
end    

lzi = lzmin:step:lzmax;
nsamp_z = length(lzi);
zi = exp(lzi);
 

laz = la*zi;
p_lz = zeros(nexp,nsamp_z);
mu_x = zeros(nexp,nsamp_z);

if ~exist('z_w'),       % Spatially invariant noise
    pg1_lz = 1./sqrt(prod(1 + laz,1));	% normalization term (depends on z, but not on Y)
    aux = exp(-0.5*V2*(1./(1+laz)));
    p_lz = aux*diag(pg1_lz);				% That gives us the conditional Gaussian density values
    										% for the observed samples and the considered samples of z
    % Compute mu_x(z) = E{x|log(z),Y}
    aux = diag(m)*(laz./(1 + laz));	% Remember: laz = la*zi
    mu_x = V.'*aux;			% Wiener estimation, for each considered sample of z
else                    % Spatially variant noise
    rep_z_w = repmat(z_w, 1, N); 
    for n_z = 1:nsamp_z,
        rep_laz = repmat(laz(:,n_z).',nexp,1);
        aux = rep_laz + rep_z_w;     % lambda*z_u + z_w
        p_lz(:,n_z) = exp(-0.5*sum(V2./aux,2))./sqrt(prod(aux,2));
        % Compute mu_x(z) = E{x|log(z),Y,z_w}
        aux = rep_laz./aux;
        mu_x(:,n_z) = (V.'.*aux)*m.';
    end
end    
                                            
                                            
[foo, ind] = max(p_lz.');	% We use ML estimation of z only for the boundaries.
clear foo
if prod(size(ind)) == 0,
	z = ones(1,size(ind,2));
else
	z = zi(ind).';				
end

clear V2 aux

% For boundary handling

uv=1+Ly;
lh=1+Lx;
dv=nblv+Ly;
rh=nblh+Lx;
ul1=ones(uv,lh);
u1=ones(uv-1,1);
l1=ones(1,lh-1);
ur1=ul1;
dl1=ul1;
dr1=ul1;
d1=u1;
r1=l1;

zM(uv:dv,lh:rh) = reshape(z,nblv,nblh);

% Propagation of the ML-estimated z to the boundaries

% a) Corners
zM(1:uv,1:lh)=zM(uv,lh)*ul1;
zM(1:uv,rh:nh)=zM(uv,rh)*ur1;
zM(dv:nv,1:lh)=zM(dv,lh)*dl1;
zM(dv:nv,rh:nh)=zM(dv,rh)*dr1;
% b) Bands
zM(1:uv-1,lh+1:rh-1)=u1*zM(uv,lh+1:rh-1);
zM(dv+1:nv,lh+1:rh-1)=d1*zM(dv,lh+1:rh-1);
zM(uv+1:dv-1,1:lh-1)=zM(uv+1:dv-1,lh)*l1;
zM(uv+1:dv-1,rh+1:nh)=zM(uv+1:dv-1,rh)*r1;

% We do scalar Wiener for the boundary coefficients
if exist('z_w'),
    x_hat = y(:,:,1).*(sx2*zM)./(sx2*zM + sig2*zW);
else        
    x_hat = y(:,:,1).*(sx2*zM)./(sx2*zM + sig2);
end


% Prior for log(z)

p_z = ones(nsamp_z,1);    % Flat log-prior (non-informative for GSM)
p_z = p_z/sum(p_z);


% Compute p(log(z)|Y) from p(Y|log(z)) and p(log(z)) (Bayes Rule)

p_lz_y = p_lz*diag(p_z);
clear p_lz
if ~optim,
   p_lz_y = (p_lz_y==max(p_lz_y')'*ones(1,size(p_lz_y,2))); 	% ML in log(z): it becomes a delta function																	% at the maximum
end    
aux = sum(p_lz_y, 2);
if any(aux==0),
    foo = aux==0;
    p_lz_y = repmat(~foo,1,nsamp_z).*p_lz_y./repmat(aux + foo,1,nsamp_z)...
        + repmat(foo,1,nsamp_z).*repmat(p_z',nexp,1); 	% Normalizing: p(log(z)|Y)
else
    p_lz_y = p_lz_y./repmat(aux,1,nsamp_z); 	% Normalizing: p(log(z)|Y)
end    
clear aux;

% Compute E{x|Y} = int_log(z){ E{x|log(z),Y} p(log(z)|Y) d(log(z)) }

aux = sum(mu_x.*p_lz_y, 2);

x_hat(1+Ly:nblv+Ly,1+Lx:nblh+Lx) = reshape(aux,nblv,nblh);

clear mu_x p_lz_y aux


function im_ext = bound_extension(im,By,Bx,type);

% im_ext = bound_extension(im,B,type);
%
% Extend an image for avoiding boundary artifacts,
%
%   By, Bx:    widths of the added stripes.
%   type:   'mirror'        Mirror extension
%           'mirror_nr':    Mirror without repeating the last pixel
%           'circular':     fft2-like
%           'zeros'

% Javier Portilla, Universidad de Granada, Jan 2004

[Ny,Nx,Nc] = size(im);

im_ext = zeros(Ny+2*By,Nx+2*Bx,Nc);
im_ext(By+1:Ny+By,Bx+1:Nx+Bx,:) = im;

if strcmp(type,'mirror'),

    im_ext(1:By,:,:) = im_ext(2*By:-1:By+1,:,:);
    im_ext(:,1:Bx,:) = im_ext(:,2*Bx:-1:Bx+1,:);
    im_ext(Ny+1+By:Ny+2*By,:,:) = im_ext(Ny+By:-1:Ny+1,:,:);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = im_ext(:,Nx+Bx:-1:Nx+1,:);
    im_ext(1:By,1:Bx,:) = im_ext(2*By:-1:By+1,2*Bx:-1:Bx+1,:);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(Ny+By:-1:Ny+1,Nx+Bx:-1:Nx+1,:);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(2*By:-1:By+1,Nx+Bx:-1:Nx+1,:);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = im_ext(Ny+By:-1:Ny+1,2*Bx:-1:Bx+1,:);

elseif strcmp(type,'mirror_nr'),    
        
    im_ext(1:By,:,:) = im_ext(2*By+1:-1:By+2,:,:);
    im_ext(:,1:Bx,:) = im_ext(:,2*Bx+1:-1:Bx+2,:);
    im_ext(Ny+1+By:Ny+2*By,:,:) = im_ext(Ny+By-1:-1:Ny,:,:);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = im_ext(:,Nx+Bx-1:-1:Nx,:);
    im_ext(1:By,1:Bx,:) = im_ext(2*By+1:-1:By+2,2*Bx+1:-1:Bx+2,:);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(Ny+By-1:-1:Ny,Nx+Bx-1:-1:Nx,:);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(2*By+1:-1:By+2,Nx+Bx-1:-1:Nx,:);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = im_ext(Ny+By-1:-1:Ny,2*Bx+1:-1:Bx+2,:);
        
elseif strcmp(type,'circular'),        
        
    im_ext(1:By,:,:) =  im_ext(Ny+1:Ny+By,:,:);
    im_ext(:,1:Bx,:) = im_ext(:,Nx+1:Nx+Bx,:);
    im_ext(Ny+1+By:Ny+2*By,:,:) = im_ext(By+1:2*By,:,:);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = im_ext(:,Bx+1:2*Bx,:);
    im_ext(1:By,1:Bx,:) = im_ext(Ny+1:Ny+By,Nx+1:Nx+Bx,:);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(By+1:2*By,Bx+1:2*Bx,:);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(Ny+1:Ny+By,Bx+1:2*Bx,:);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = im_ext(By+1:2*By,Nx+1:Nx+Bx,:);
   
end    
        

% M = MEAN2(MTX)
%
% Sample mean of a matrix.

function res = mean2(mtx)

res = mean(mean(mtx));

function [pyr,pind] = buildWUpyr(im, Nsc, daub_order);

% [PYR, INDICES] = buildWUpyr(IM, HEIGHT, DAUB_ORDER)
% 
% Construct a separable undecimated orthonormal QMF/wavelet pyramid
% on matrix (or vector) IM.
% 
% HEIGHT specifies the number of pyramid levels to build. Default
% is maxPyrHt(IM,FILT).  You can also specify 'auto' to use this value.
% 
% DAUB_ORDER: specifies the order of the daubechies wavelet filter used
% 
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband. 

% JPM, Univ. de Granada, 03/2003, based on Rice Wavelet Toolbox 
% function "mrdwt" and on Matlab Pyrtools from Eero Simoncelli.

if Nsc < 1,
    display('Error: Number of scales must be >=1.');
else,   

Nor = 3; % fixed number of orientations;
h = daubcqf(daub_order);
[lpr,yh] = mrdwt(im, h, Nsc+1); % performs the decomposition

[Ny,Nx] = size(im);

% Reorganize the output, forcing the same format as with buildFullSFpyr2

pyr = [];
pind = zeros((Nsc+1)*Nor+2,2);    % Room for a "virtual" high pass residual, for compatibility
nband = 1;
for nsc = 1:Nsc+1,
    for nor = 1:Nor,
        nband = nband + 1;
        band = yh(:,(nband-2)*Nx+1:(nband-1)*Nx);
        sh = (daub_order/2 - 1)*2^nsc;  % approximate phase compensation
        if nor == 1,        % horizontal
            band = shift(band, [sh 2^(nsc-1)]);
        elseif nor == 2,    % vertical
            band = shift(band, [2^(nsc-1) sh]);
        else
            band = shift(band, [sh sh]);    % diagonal
        end    
        if nsc>2,
            band = real(shrink(band,2^(nsc-2)));  % The low freq. bands are shrunk in the freq. domain
        end
        pyr = [pyr; vector(band)];
        pind(nband,:) = size(band);
    end    
end            

band = lpr;
band = shrink(band,2^Nsc);
pyr = [pyr; vector(band)];
pind(nband+1,:) = size(band);


end


function [h_0,h_1] = daubcqf(N,TYPE)
%    [h_0,h_1] = daubcqf(N,TYPE); 
%
%    Function computes the Daubechies' scaling and wavelet filters
%    (normalized to sqrt(2)).
%
%    Input: 
%       N    : Length of filter (must be even)
%       TYPE : Optional parameter that distinguishes the minimum phase,
%              maximum phase and mid-phase solutions ('min', 'max', or
%              'mid'). If no argument is specified, the minimum phase
%              solution is used.
%
%    Output: 
%       h_0 : Minimal phase Daubechies' scaling filter 
%       h_1 : Minimal phase Daubechies' wavelet filter 
%
%    Example:
%       N = 4;
%       TYPE = 'min';
%       [h_0,h_1] = daubcqf(N,TYPE)
%       h_0 = 0.4830 0.8365 0.2241 -0.1294
%       h_1 = 0.1294 0.2241 -0.8365 0.4830
%
%    Reference: "Orthonormal Bases of Compactly Supported Wavelets",
%                CPAM, Oct.89 
%

%File Name: daubcqf.m
%Last Modification Date: 01/02/96	15:12:57
%Current Version: daubcqf.m	2.4
%File Creation Date: 10/10/88
%Author: Ramesh Gopinath  <ramesh@dsp.rice.edu>
%
%Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
%Created by Ramesh Gopinath, Department of ECE, Rice University. 
%
%This software is distributed and licensed to you on a non-exclusive 
%basis, free-of-charge. Redistribution and use in source and binary forms, 
%with or without modification, are permitted provided that the following 
%conditions are met:
%
%1. Redistribution of source code must retain the above copyright notice, 
%   this list of conditions and the following disclaimer.
%2. Redistribution in binary form must reproduce the above copyright notice, 
%   this list of conditions and the following disclaimer in the 
%   documentation and/or other materials provided with the distribution.
%3. All advertising materials mentioning features or use of this software 
%   must display the following acknowledgment: This product includes 
%   software developed by Rice University, Houston, Texas and its contributors.
%4. Neither the name of the University nor the names of its contributors 
%   may be used to endorse or promote products derived from this software 
%   without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY WILLIAM MARSH RICE UNIVERSITY, HOUSTON, TEXAS, 
%AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
%BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
%FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE UNIVERSITY 
%OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
%OR BUSINESS INTERRUPTIONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
%WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
%OTHERWISE), PRODUCT LIABILITY, OR OTHERWISE ARISING IN ANY WAY OUT OF THE 
%USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%For information on commercial licenses, contact Rice University's Office of 
%Technology Transfer at techtran@rice.edu or (713) 348-6173

if(nargin < 2),
  TYPE = 'min';
end;
if(rem(N,2) ~= 0),
  error('No Daubechies filter exists for ODD length');
end;
K = N/2;
a = 1;
p = 1;
q = 1;
h_0 = [1 1];
for j  = 1:K-1,
  a = -a * 0.25 * (j + K - 1)/j;
  h_0 = [0 h_0] + [h_0 0];
  p = [0 -p] + [p 0];
  p = [0 -p] + [p 0];
  q = [0 q 0] + a*p;
end;
q = sort(roots(q));
qt = q(1:K-1);
if TYPE=='mid',
  if rem(K,2)==1,  
    qt = q([1:4:N-2 2:4:N-2]);
  else
    qt = q([1 4:4:K-1 5:4:K-1 N-3:-4:K N-4:-4:K]);
  end;
end;
h_0 = conv(h_0,real(poly(qt)));
h_0 = sqrt(2)*h_0/sum(h_0); 	%Normalize to sqrt(2);
if(TYPE=='max'),
  h_0 = fliplr(h_0);
end;
if(abs(sum(h_0 .^ 2))-1 > 1e-4) 
  error('Numerically unstable for this value of "N".');
end;
h_1 = rot90(h_0,2);
h_1(1:2:N)=-h_1(1:2:N);


% [RES] = shift(MTX, OFFSET)
% 
% Circular shift 2D matrix samples by OFFSET (a [Y,X] 2-vector),
% such that  RES(POS) = MTX(POS-OFFSET).

function res = shift(mtx, offset)

dims = size(mtx);

offset = mod(-offset,dims);

res = [ mtx(offset(1)+1:dims(1), offset(2)+1:dims(2)),  ...
          mtx(offset(1)+1:dims(1), 1:offset(2));        ...
        mtx(1:offset(1), offset(2)+1:dims(2)),          ...
	  mtx(1:offset(1), 1:offset(2)) ];

  
% [VEC] = vector(MTX)
% 
% Pack elements of MTX into a column vector.  Same as VEC = MTX(:)
% Previously named "vectorize" (changed to avoid overlap with Matlab's
% "vectorize" function).

function vec = vector(mtx)

vec = mtx(:);


function ts = shrink(t,f)

%	im_shr = shrink(im0, f)
%
% It shrinks (spatially) an image into a factor f
% in each dimension. It does it by cropping the
% Fourier transform of the image.

% JPM, 5/1/95.

% Revised so it can work also with exponents of 3 factors: JPM 5/2003


[my,mx] = size(t);
T = fftshift(fft2(t))/f^2;
Ts = zeros(my/f,mx/f);

cy = ceil(my/2);
cx = ceil(mx/2);
evenmy = (my/2==floor(my/2));
evenmx = (mx/2==floor(mx/2));

y1 = cy + 2*evenmy - floor(my/(2*f));
y2 = cy + floor(my/(2*f));
x1 = cx + 2*evenmx - floor(mx/(2*f));
x2 = cx + floor(mx/(2*f));

Ts(1+evenmy:my/f,1+evenmx:mx/f)=T(y1:y2,x1:x2);
if evenmy,
    Ts(1+evenmy:my/f,1)=(T(y1:y2,x1-1)+T(y1:y2,x2+1))/2;
end
if evenmx,
    Ts(1,1+evenmx:mx/f)=(T(y1-1,x1:x2)+T(y2+1,x1:x2))/2;
end
if evenmy & evenmx,
    Ts(1,1)=(T(y1-1,x1-1)+T(y1-1,x2+1)+T(y2+1,x1-1)+T(y2+1,x2+1))/4;
end    
Ts = fftshift(Ts);
Ts = shift(Ts, [1 1] - [evenmy evenmx]);
ts = ifft2(Ts);

% RES = pyrBand(PYR, INDICES, BAND_NUM)
%
% Access a subband from a pyramid (gaussian, laplacian, QMF/wavelet, 
% or steerable).  Subbands are numbered consecutively, from finest
% (highest spatial frequency) to coarsest (lowest spatial frequency).

% Eero Simoncelli, 6/96.

function res =  pyrBand(pyr, pind, band)

res = reshape( pyr(pyrBandIndices(pind,band)), pind(band,1), pind(band,2) );


% RES = pyrBandIndices(INDICES, BAND_NUM)
%
% Return indices for accessing a subband from a pyramid 
% (gaussian, laplacian, QMF/wavelet, steerable).

% Eero Simoncelli, 6/96.

function indices =  pyrBandIndices(pind,band)

if ((band > size(pind,1)) | (band < 1))
  error(sprintf('BAND_NUM must be between 1 and number of pyramid bands (%d).', ...
      size(pind,1)));
end

if (size(pind,2) ~= 2)
  error('INDICES must be an Nx2 matrix indicating the size of the pyramid subbands');
end

ind = 1;
for l=1:band-1
  ind = ind + prod(pind(l,:));
end

indices = ind:ind+prod(pind(band,:))-1;


function res = reconWUpyr(pyr, pind, daub_order);

% RES = reconWUpyr(PYR, INDICES, DAUB_ORDER)
 
% Reconstruct image from its separable undecimated orthonormal QMF/wavelet pyramid
% representation, as created by buildWUpyr.
%
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  
% 
% DAUB_ORDER: specifies the order of the daubechies wavelet filter used
 
% JPM, Univ. de Granada, 03/2003, based on Rice Wavelet Toolbox 
% functions "mrdwt" and  "mirdwt", and on Matlab Pyrtools from Eero Simoncelli.


Nor = 3;
Nsc = (size(pind,1)-2)/Nor-1;
h = daubcqf(daub_order);

yh = [];

nband = 1;
last = prod(pind(1,:)); % empty "high pass residual band" for compatibility with full steerpyr 2
for nsc = 1:Nsc+1,  % The number of scales corresponds to the number of pyramid levels (also for compatibility)
    for nor = 1:Nor,
        nband = nband +1;
        first = last + 1;
        last = first + prod(pind(nband,:)) - 1;
        band = pyrBand(pyr,pind,nband);
        sh = (daub_order/2 - 1)*2^nsc;  % approximate phase compensation
        if nsc > 2,
            band = expand(band, 2^(nsc-2));
        end   
        if nor == 1,        % horizontal
            band = shift(band, [-sh -2^(nsc-1)]);
        elseif nor == 2,    % vertical
            band = shift(band, [-2^(nsc-1) -sh]);
        else
            band = shift(band, [-sh -sh]);    % diagonal
        end    
        yh = [yh band];    
    end    
end    

nband = nband + 1;
band = pyrBand(pyr,pind,nband);
lpr = expand(band,2^Nsc);

res= mirdwt(lpr,yh,h,Nsc+1);

function te = expand(t,f)

%	im_exp = expand(im0, f)
%
% It expands (spatially) an image into a factor f
% in each dimension. It does it filling in with zeros
% the expanded Fourier domain.

% JPM, 5/1/95.

% Revised so it can work also with exponents of 3 factors: JPM 5/2003

[my mx] = size(t);
my = f*my;
mx = f*mx;
Te = zeros(my,mx);
T = f^2*fftshift(fft2(t));

cy = ceil(my/2);
cx = ceil(mx/2);
evenmy = (my/2==floor(my/2));
evenmx = (mx/2==floor(mx/2));

y1 = cy + 2*evenmy - floor(my/(2*f));
y2 = cy + floor(my/(2*f));
x1 = cx + 2*evenmx - floor(mx/(2*f));
x2 = cx + floor(mx/(2*f));

Te(y1:y2,x1:x2)=T(1+evenmy:my/f,1+evenmx:mx/f);
if evenmy,
    Te(y1-1,x1:x2)=T(1,2:mx/f)/2;
    Te(y2+1,x1:x2)=((T(1,mx/f:-1:2)/2)').';
end
if evenmx,
    Te(y1:y2,x1-1)=T(2:my/f,1)/2;
    Te(y1:y2,x2+1)=((T(my/f:-1:2,1)/2)').';
end
if evenmx & evenmy,
    esq=T(1,1)/4;
    Te(y1-1,x1-1)=esq;
    Te(y1-1,x2+1)=esq;
    Te(y2+1,x1-1)=esq;
    Te(y2+1,x2+1)=esq;
end    
Te=fftshift(Te);
Te = shift(Te, [1 1] - [evenmy evenmx]);
te=ifft2(Te);


% RES = innerProd(MTX)
%
% Compute (MTX' * MTX) efficiently (i.e., without copying the matrix)

function res = innerProd(mtx)

%% NOTE: THIS CODE SHOULD NOT BE USED! (MEX FILE IS CALLED INSTEAD)

% fprintf(1,'WARNING: You should compile the MEX version of "innerProd.c",\n         found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster.\n');

res = mtx' * mtx;
