function y = perform_pyramid_transform_simoncelli(x, Jmin,options)

% perform_pyramid_transform_simoncelli - Laplacian pyramidal transform
%
%   y = perform_pyramid_transform_simoncelli(x, Jmin);
%
%   This is just a convenient wrapper to the original steerable 
%   matlab toolbox of Simoncelli that can be downloaded from
%       http://www.cns.nyu.edu/~eero/STEERPYR/
%   If you want a fast implementation of this function, you
%   have to install Simoncelli toolbox (to use mex files).
%
%   It provide a simpler interface that directly output a cell
%   array of images. Usage :
%
%   M = load_image('lena');
%   MS = perform_pyramid_transform_simoncelli(M, 3); % synthesis
%   M1 = perform_pyramid_transform_simoncelli(MS);   % reconstruction
%   
%   Copyright (c) 2005 Gabriel Peyré

if nargin<3
    options.null = 0;
end

filts = 'binom5';

if ~iscell(x)
    if nargin<2
        Jmin = 4;
    end
    Jmax = log2(size(x,1))-1;
    nbr_bands = Jmax - Jmin + 2;
    % fwd transform
    [pyr,pind] = buildLpyr(x, nbr_bands, filts, filts, 'reflect1');
    % copy into cell array    
    y = {};
    for k=1:size(pind, 1)
        indices =  pyrBandIndices(pind,k);
        L = length(indices);    % length of this scale
        y{k} = reshape( pyr(indices), sqrt(L), sqrt(L) );
    end
else
    n = size(x{1},1);
    Jmax = log2(n)-1;
    nbr_bands = length(x);
    % copy from cell array
    pind = n ./ 2.^(0:nbr_bands-1);
    pind = [pind(:), pind(:)];
    % build the matrix
    n = sum( prod(pind,2) );
    pyr = zeros(n, 1);
    for k=1:size(pind, 1)
        indices =  pyrBandIndices(pind,k);
        L = length(indices);    % length of this scale
        pyr(indices) = x{k}(:);
    end
    % bwd transform
    y = reconLpyr(pyr, pind, 'all', filts); 
end





% [PYR, INDICES] = buildLpyr(IM, HEIGHT, FILT1, FILT2, EDGES)
%
% Construct a Laplacian pyramid on matrix (or vector) IM.
%
% HEIGHT (optional) specifies the number of pyramid levels to build. Default
% is 1+maxPyrHt(size(IM),size(FILT)).  You can also specify 'auto' to
% use this value.
%
% FILT1 (optional) can be a string naming a standard filter (see
% namedFilter), or a vector which will be used for (separable)
% convolution.  Default = 'binom5'.  FILT2 specifies the "expansion"
% filter (default = filt1).  EDGES specifies edge-handling, and
% defaults to 'reflect1' (see corrDn).
%
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  This is compatible with the MatLab Wavelet toolbox.

% Eero Simoncelli, 6/96.

function [pyr,pind] = buildLpyr(im, ht, filt1, filt2, edges)

if (nargin < 1)
  error('First argument (IM) is required');
end

im_sz = size(im);

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('filt1') ~= 1)
  filt1 = 'binom5';
end
 
if isstr(filt1)
  filt1 = namedFilter(filt1);
end

if ( (size(filt1,1) > 1) & (size(filt1,2) > 1) )
  error('FILT1 should be a 1D filter (i.e., a vector)');
else
  filt1 = filt1(:);
end

if (exist('filt2') ~= 1)
  filt2 = filt1;
end

if isstr(filt2)
  filt2 = namedFilter(filt2);
end

if ( (size(filt2,1) > 1) & (size(filt2,2) > 1) )
  error('FILT2 should be a 1D filter (i.e., a vector)');
else
  filt2 = filt2(:);
end

max_ht = 1 + maxPyrHt(im_sz, max(size(filt1,1), size(filt2,1)));
if ( (exist('ht') ~= 1) | (ht == 'auto') )
  ht = max_ht;
else
  if (ht > max_ht)
    error(sprintf('Cannot build pyramid higher than %d levels.',max_ht));
  end
end

if (exist('edges') ~= 1)
  edges= 'reflect1';
end

%------------------------------------------------------------

if (ht <= 1)

  pyr = im(:);
  pind = im_sz;

else

  if (im_sz(2) == 1)
    lo2 = corrDn(im, filt1, edges, [2 1], [1 1]);
  elseif (im_sz(1) == 1)
    lo2 = corrDn(im, filt1', edges, [1 2], [1 1]);
  else
    lo = corrDn(im, filt1', edges, [1 2], [1 1]);
    int_sz = size(lo);
    lo2 = corrDn(lo, filt1, edges, [2 1], [1 1]);
  end

  [npyr,nind] = buildLpyr(lo2, ht-1, filt1, filt2, edges);

  if (im_sz(1) == 1)
    hi2 = upConv(lo2, filt2', edges, [1 2], [1 1], im_sz);
  elseif (im_sz(2) == 1)
    hi2 = upConv(lo2, filt2, edges, [2 1], [1 1], im_sz);
  else
    hi = upConv(lo2, filt2, edges, [2 1], [1 1], int_sz);
    hi2 = upConv(hi, filt2', edges, [1 2], [1 1], im_sz);
  end

  hi2 = im - hi2;

  pyr = [hi2(:); npyr];
  pind = [im_sz; nind];

end
  


% KERNEL = NAMED_FILTER(NAME)
%
% Some standard 1D filter kernels.  These are scaled such that
% their L2-norm is 1.0.
%
%  binomN             - binomial coefficient filter of order N-1
%  haar:              - Haar wavelet.
%  qmf8, qmf12, qmf16 - Symmetric Quadrature Mirror Filters [Johnston80]
%  daub2,daub3,daub4  - Daubechies wavelet [Daubechies88].
%  qmf5, qmf9, qmf13: - Symmetric Quadrature Mirror Filters [Simoncelli88,Simoncelli90]
%
%  See bottom of file for full citations.

% Eero Simoncelli, 6/96.

function [kernel] = namedFilter(name)

if strcmp(name(1:min(5,size(name,2))), 'binom')
  kernel = sqrt(2) * binomialFilter(str2num(name(6:size(name,2))));
elseif strcmp(name,'qmf5')
  kernel = [-0.076103 0.3535534 0.8593118 0.3535534 -0.076103]';
elseif strcmp(name,'qmf9')
  kernel = [0.02807382 -0.060944743 -0.073386624 0.41472545 0.7973934 ...
      0.41472545 -0.073386624 -0.060944743 0.02807382]';
elseif strcmp(name,'qmf13')
  kernel = [-0.014556438 0.021651438 0.039045125 -0.09800052 ...
	  -0.057827797 0.42995453 0.7737113 0.42995453 -0.057827797 ...
	  -0.09800052 0.039045125 0.021651438 -0.014556438]';
elseif strcmp(name,'qmf8')
  kernel = sqrt(2) * [0.00938715 -0.07065183 0.06942827 0.4899808 ...
    0.4899808 0.06942827 -0.07065183 0.00938715 ]';
elseif strcmp(name,'qmf12')
  kernel = sqrt(2) * [-0.003809699 0.01885659 -0.002710326 -0.08469594 ...
	0.08846992 0.4843894 0.4843894 0.08846992 -0.08469594 -0.002710326 ...
	0.01885659 -0.003809699 ]';
elseif strcmp(name,'qmf16')
  kernel = sqrt(2) * [0.001050167 -0.005054526 -0.002589756 0.0276414 -0.009666376 ...
	-0.09039223 0.09779817 0.4810284 0.4810284 0.09779817 -0.09039223 -0.009666376 ...
	0.0276414 -0.002589756 -0.005054526 0.001050167 ]';
elseif strcmp(name,'haar')
  kernel = [1 1]' / sqrt(2);
elseif strcmp(name,'daub2')
  kernel = [0.482962913145 0.836516303738 0.224143868042 -0.129409522551]';
elseif strcmp(name,'daub3')
  kernel = [0.332670552950 0.806891509311 0.459877502118 -0.135011020010 ...
	-0.085441273882  0.035226291882]';
elseif strcmp(name,'daub4')
  kernel = [0.230377813309 0.714846570553 0.630880767930 -0.027983769417 ...
	-0.187034811719 0.030841381836 0.032883011667 -0.010597401785]';
elseif strcmp(name,'gauss5')  % for backward-compatibility
  kernel = sqrt(2) * [0.0625 0.25 0.375 0.25 0.0625]';
elseif strcmp(name,'gauss3')  % for backward-compatibility
  kernel = sqrt(2) * [0.25 0.5 0.25]';
else
  error(sprintf('Bad filter name: %s\n',name));
end
  
% [Johnston80] - J D Johnston, "A filter family designed for use in quadrature 
%    mirror filter banks", Proc. ICASSP, pp 291-294, 1980.
%
% [Daubechies88] - I Daubechies, "Orthonormal bases of compactly supported wavelets",
%    Commun. Pure Appl. Math, vol. 42, pp 909-996, 1988.
%
% [Simoncelli88] - E P Simoncelli,  "Orthogonal sub-band image transforms",
%     PhD Thesis, MIT Dept. of Elec. Eng. and Comp. Sci. May 1988.
%     Also available as: MIT Media Laboratory Vision and Modeling Technical 
%     Report #100.
%
% [Simoncelli90] -  E P Simoncelli and E H Adelson, "Subband image coding",
%    Subband Transforms, chapter 4, ed. John W Woods, Kluwer Academic 
%    Publishers,  Norwell, MA, 1990, pp 143--192.


% KERNEL = binomialFilter(size)
%
% Returns a vector of binomial coefficients of order (size-1) .

% Eero Simoncelli, 2/97.

function [kernel] = binomialFilter(sz)

if (sz < 2)
  error('size argument must be larger than 1');
end

kernel = [0.5 0.5]';

for n=1:sz-2
  kernel = conv([0.5 0.5]', kernel);
end
  
% HEIGHT = maxPyrHt(IMSIZE, FILTSIZE)
%
% Compute maximum pyramid height for given image and filter sizes.
% Specifically: the number of corrDn operations that can be sequentially
% performed when subsampling by a factor of 2.

% Eero Simoncelli, 6/96.

function height = maxPyrHt(imsz, filtsz)

imsz = imsz(:);
filtsz = filtsz(:);

if any(imsz == 1) % 1D image
  imsz = prod(imsz);
  filtsz = prod(filtsz);
elseif any(filtsz == 1)              % 2D image, 1D filter
  filtsz = [filtsz(1); filtsz(1)];
end

if any(imsz < filtsz)
  height = 0;
else
  height = 1 + maxPyrHt( floor(imsz/2), filtsz ); 
end


% RES = corrDn(IM, FILT, EDGES, STEP, START, STOP)
%
% Compute correlation of matrices IM with FILT, followed by
% downsampling.  These arguments should be 1D or 2D matrices, and IM
% must be larger (in both dimensions) than FILT.  The origin of filt
% is assumed to be floor(size(filt)/2)+1.
% 
% EDGES is a string determining boundary handling:
%    'circular' - Circular convolution
%    'reflect1' - Reflect about the edge pixels
%    'reflect2' - Reflect, doubling the edge pixels
%    'repeat'   - Repeat the edge pixels
%    'zero'     - Assume values of zero outside image boundary
%    'extend'   - Reflect and invert
%    'dont-compute' - Zero output when filter overhangs input boundaries
%
% Downsampling factors are determined by STEP (optional, default=[1 1]), 
% which should be a 2-vector [y,x].
% 
% The window over which the convolution occurs is specfied by START 
% (optional, default=[1,1], and STOP (optional, default=size(IM)).
% 
% NOTE: this operation corresponds to multiplication of a signal
% vector by a matrix whose rows contain copies of the FILT shifted by
% multiples of STEP.  See upConv.m for the operation corresponding to
% the transpose of this matrix.

% Eero Simoncelli, 6/96, revised 2/97.

function res = corrDn(im, filt, edges, step, start, stop)

%% NOTE: THIS CODE IS NOT ACTUALLY USED! (MEX FILE IS CALLED INSTEAD)

% fprintf(1,'WARNING: You should compile the MEX version of "corrDn.c",\n         found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster, and provides more boundary-handling options.\n');

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('edges') == 1) 
  if (strcmp(edges,'reflect1') ~= 1)
    warning('Using REFLECT1 edge-handling (use MEX code for other options).');
  end
end

if (exist('step') ~= 1)
	step = [1,1];
end	

if (exist('start') ~= 1)
	start = [1,1];
end	

if (exist('stop') ~= 1)
	stop = size(im);
end	

%------------------------------------------------------------

% Reverse order of taps in filt, to do correlation instead of convolution
filt = filt(size(filt,1):-1:1,size(filt,2):-1:1);

tmp = rconv2(im,filt);
res = tmp(start(1):step(1):stop(1),start(2):step(2):stop(2));


% RES = RCONV2(MTX1, MTX2, CTR)
%
% Convolution of two matrices, with boundaries handled via reflection
% about the edge pixels.  Result will be of size of LARGER matrix.
% 
% The origin of the smaller matrix is assumed to be its center.
% For even dimensions, the origin is determined by the CTR (optional) 
% argument:
%      CTR   origin
%       0     DIM/2      (default)
%       1     (DIM/2)+1  

% Eero Simoncelli, 6/96.

function c = rconv2(a,b,ctr)

if (exist('ctr') ~= 1)
  ctr = 0;
end

if (( size(a,1) >= size(b,1) ) & ( size(a,2) >= size(b,2) ))
    large = a; small = b;
elseif  (( size(a,1) <= size(b,1) ) & ( size(a,2) <= size(b,2) ))
    large = b; small = a;
else
  error('one arg must be larger than the other in both dimensions!');
end

ly = size(large,1);
lx = size(large,2);
sy = size(small,1);
sx = size(small,2);

%% These values are one less than the index of the small mtx that falls on 
%% the border pixel of the large matrix when computing the first 
%% convolution response sample:
sy2 = floor((sy+ctr-1)/2);
sx2 = floor((sx+ctr-1)/2);

% pad with reflected copies
clarge = [ 
    large(sy-sy2:-1:2,sx-sx2:-1:2), large(sy-sy2:-1:2,:), ...
	large(sy-sy2:-1:2,lx-1:-1:lx-sx2); ...
    large(:,sx-sx2:-1:2),    large,   large(:,lx-1:-1:lx-sx2); ...
    large(ly-1:-1:ly-sy2,sx-sx2:-1:2), ...
      large(ly-1:-1:ly-sy2,:), ...
      large(ly-1:-1:ly-sy2,lx-1:-1:lx-sx2) ];

c = conv2(clarge,small,'valid');

% RES = upConv(IM, FILT, EDGES, STEP, START, STOP, RES)
%
% Upsample matrix IM, followed by convolution with matrix FILT.  These
% arguments should be 1D or 2D matrices, and IM must be larger (in
% both dimensions) than FILT.  The origin of filt
% is assumed to be floor(size(filt)/2)+1.
%
% EDGES is a string determining boundary handling:
%    'circular' - Circular convolution
%    'reflect1' - Reflect about the edge pixels
%    'reflect2' - Reflect, doubling the edge pixels
%    'repeat'   - Repeat the edge pixels
%    'zero'     - Assume values of zero outside image boundary
%    'extend'   - Reflect and invert
%    'dont-compute' - Zero output when filter overhangs OUTPUT boundaries
%
% Upsampling factors are determined by STEP (optional, default=[1 1]),
% a 2-vector [y,x].
% 
% The window over which the convolution occurs is specfied by START 
% (optional, default=[1,1], and STOP (optional, default = 
% step .* (size(IM) + floor((start-1)./step))).
%
% RES is an optional result matrix.  The convolution result will be 
% destructively added into this matrix.  If this argument is passed, the 
% result matrix will not be returned. DO NOT USE THIS ARGUMENT IF 
% YOU DO NOT UNDERSTAND WHAT THIS MEANS!!
% 
% NOTE: this operation corresponds to multiplication of a signal
% vector by a matrix whose columns contain copies of the time-reversed
% (or space-reversed) FILT shifted by multiples of STEP.  See corrDn.m
% for the operation corresponding to the transpose of this matrix.

% Eero Simoncelli, 6/96.  revised 2/97.

function result = upConv(im,filt,edges,step,start,stop,res)

%% THIS CODE IS NOT ACTUALLY USED! (MEX FILE IS CALLED INSTEAD)

% fprintf(1,'WARNING: You should compile the MEX version of "upConv.c",\n         found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster, and provides more boundary-handling options.\n');

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('edges') == 1) 
  if (strcmp(edges,'reflect1') ~= 1)
    warning('Using REFLECT1 edge-handling (use MEX code for other options).');
  end
end

if (exist('step') ~= 1)
  step = [1,1];
end	

if (exist('start') ~= 1)
  start = [1,1];
end	

% A multiple of step
if (exist('stop') ~= 1)
  stop = step .* (floor((start-ones(size(start)))./step)+size(im))
end	

if ( ceil((stop(1)+1-start(1)) / step(1)) ~= size(im,1) )
  error('Bad Y result dimension');
end
if ( ceil((stop(2)+1-start(2)) / step(2)) ~= size(im,2) )
  error('Bad X result dimension');
end

if (exist('res') ~= 1)
  res = zeros(stop-start+1);
end	

%------------------------------------------------------------

tmp = zeros(size(res));
tmp(start(1):step(1):stop(1),start(2):step(2):stop(2)) = im;

result = rconv2(tmp,filt) + res;


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


% RES = reconLpyr(PYR, INDICES, LEVS, FILT2, EDGES)
%
% Reconstruct image from Laplacian pyramid, as created by buildLpyr.
%
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  This is compatible with the MatLab Wavelet toolbox.
%
% LEVS (optional) should be a list of levels to include, or the string
% 'all' (default).  The finest scale is number 1.  The lowpass band
% corresponds to lpyrHt(INDICES)+1.
%
% FILT2 (optional) can be a string naming a standard filter (see
% namedFilter), or a vector which will be used for (separable)
% convolution.  Default = 'binom5'.  EDGES specifies edge-handling,
% and defaults to 'reflect1' (see corrDn).

% Eero Simoncelli, 6/96

function res = reconLpyr(pyr, ind, levs, filt2, edges)

if (nargin < 2)
  error('First two arguments (PYR, INDICES) are required');
end
  
%%------------------------------------------------------------
%% DEFAULTS:

if (exist('levs') ~= 1)
  levs = 'all';
end

if (exist('filt2') ~= 1)
  filt2 = 'binom5';
end

if (exist('edges') ~= 1)
  edges= 'reflect1';
end
%%------------------------------------------------------------

maxLev =  1+lpyrHt(ind);
if strcmp(levs,'all')
  levs = [1:maxLev]';
else
  if (any(levs > maxLev))
    error(sprintf('Level numbers must be in the range [1, %d].', maxLev));
  end
  levs = levs(:);
end

if isstr(filt2)
  filt2 = namedFilter(filt2);
end

filt2 = filt2(:);
res_sz = ind(1,:);

if any(levs > 1)

  int_sz = [ind(1,1), ind(2,2)];
  
  nres = reconLpyr( pyr(prod(res_sz)+1:size(pyr,1)), ...
      ind(2:size(ind,1),:), levs-1, filt2, edges);
  
  if (res_sz(1) == 1)
    res = upConv(nres, filt2', edges, [1 2], [1 1], res_sz);
  elseif (res_sz(2) == 1)
    res = upConv(nres, filt2, edges, [2 1], [1 1], res_sz);
  else
    hi = upConv(nres, filt2, edges, [2 1], [1 1], int_sz);
    res = upConv(hi, filt2', edges, [1 2], [1 1], res_sz);
  end

else
  
  res = zeros(res_sz);

end

if any(levs == 1)
  res = res + pyrBand(pyr,ind,1);
end


% [HEIGHT] = lpyrHt(INDICES)
%
% Compute height of Laplacian pyramid with given its INDICES matrix.
% See buildLpyr.m

% Eero Simoncelli, 6/96.

function [ht] =  lpyrHt(pind)

% Don't count lowpass residual band
ht = size(pind,1)-1;


% RES = pyrBand(PYR, INDICES, BAND_NUM)
%
% Access a subband from a pyramid (gaussian, laplacian, QMF/wavelet, 
% or steerable).  Subbands are numbered consecutively, from finest
% (highest spatial frequency) to coarsest (lowest spatial frequency).

% Eero Simoncelli, 6/96.

function res =  pyrBand(pyr, pind, band)

res = reshape( pyr(pyrBandIndices(pind,band)), pind(band,1), pind(band,2) );
