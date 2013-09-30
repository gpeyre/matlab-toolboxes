function x = perform_histogram_matching(x, y, options)

% perform_histogram_matching - match the histogram of two image.
%
%   x = perform_histogram_matching(x, y, nb_bins);
% or
%   x = perform_histogram_matching(x, y, options);
%
%   Perform an equalization of x so that it histogram 
%   matches the histogram of y.
%
%   Works also for vector valued (ie 3D matrix) images.
%   For color images, ie (n,p,3) sized image, the matching is performed
%       *  in RVB space if options.match_ycbcr=0
%       *  in YcBCr space if options.match_ycbcr=1 (default, works better)
%
%   Copyright (c) 2005 Gabriel Peyre

if nargin>2 && not(isstruct(options))
    opt = options;
    clear options;
    options.nb_bins = opt;
end

options.null = 0;
nb_bins = getoptions(options, 'nb_bins', 100);
match_ycbcr = getoptions(options, 'match_ycbcr', 1);
absval = getoptions(options, 'absval', 0);
cols = getoptions(options, 'cols', 0);
rows = getoptions(options, 'rows', 0);

if cols && rows
    error('You cannote specify both cols and rows');
end

if cols && size(x,2)>1
    if not(size(x,2)==size(y,2))
        error('x and y must have same number of columns');
    end
    if size(x,3)>1 || size(y,3)>1
        error('options.cols does not works for color images');
    end
    for i=1:size(x,2)
        x(:,i) = perform_histogram_matching(x(:,i),y(:,i),options);
    end
    return;
end
if cols && size(x,1)>1
    if not(size(x,1)==size(y,1))
        error('x and y must have same number of rows');
    end
    if size(x,3)>1 || size(y,3)>1
        error('options.rows does not works for color images');
    end
    for i=1:size(x,1)
        x(i,:) = perform_histogram_matching(x(i,:),y(i,:),options);
    end
    return;
end


if iscell(x)
   if ~(length(x)==length(y))
       error('For cell arrays, src and tgt must have the same length.');
   end
   for i=1:length(x)
       if length(find(isinf(x{i})))==0
            x{i} = perform_histogram_matching(x{i}, y{i}, options);
       else
			x{i} = x{i};
       end
   end
   return;
end


if size(x,3)>1
    if size(x,3)~=size(y,3)
        error('Src and tgt images must have the same number of components.');
    end
    
    if size(x,3)==3 && match_ycbcr
        x = rgb2ycbcr(x);
        y = rgb2ycbcr(y);
    end
    % match each color
    for i=1:size(y,3)
        x(:,:,i) = perform_histogram_matching(x(:,:,i), y(:,:,i), options);
    end
    if size(x,3)==3 && match_ycbcr
        x = ycbcr2rgb(x);
    end
    return;
end


if sum(abs(y(:) - mean(y(:))))<1e-8
    % histo is crashing for constant signals
    x = x*0 + mean(y(:));
    return;
end

sx = size(x);
x = x(:);
y = y(:);

if absval
    s = sign(x);
    x = abs(x);
    y = abs(y);
end


% compute transformed histograms
% [N,X] = histo(y, nb_bins);
[N,X] = histo_slow(y, nb_bins);

x = histoMatch(x, N, X);

I = find(isnan(x(:)));
x(I) = max(y(:));

if absval
    x = x .* s;
end


x = reshape(x,sx);


% [N,X] = histo_slow(MTX, nbinsOrBinsize, binCenter);
%
% Compute a histogram of (all) elements of MTX.  N contains the histogram
% counts, X is a vector containg the centers of the histogram bins.
%
% nbinsOrBinsize (optional, default = 101) specifies either
% the number of histogram bins, or the negative of the binsize.
%
% binCenter (optional, default = mean2(MTX)) specifies a center position
% for (any one of) the histogram bins.
%
% How does this differ from MatLab's HIST function?  This function:
%   - allows uniformly spaced bins only.
%   +/- operates on all elements of MTX, instead of columnwise.
%   + is much faster (approximately a factor of 80 on my machine).
%   + allows specification of number of bins OR binsize.  Default=101 bins.
%   + allows (optional) specification of binCenter.

% Eero Simoncelli, 3/97.

function [N, X] = histo_slow(mtx, nbins, binCtr)

%% NOTE: THIS CODE IS NOT ACTUALLY USED! (MEX FILE IS CALLED INSTEAD)

% fprintf(1,'WARNING: You should compile the MEX version of "histo.c",\n         found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster.\n');

mtx = mtx(:);

%------------------------------------------------------------
%% OPTIONAL ARGS:

[mn,mx] = range2(mtx);

if (exist('binCtr') ~= 1) 
  binCtr =  mean(mtx);
end

if (exist('nbins') == 1) 
  if (nbins < 0)
    binSize = -nbins;
  else
    binSize = ((mx-mn)/nbins);
    tmpNbins = round((mx-binCtr)/binSize) - round((mn-binCtr)/binSize);
    if (tmpNbins ~= nbins)
      warning('Using %d bins instead of requested number (%d)',tmpNbins,nbins);
    end
  end
else
  binSize = ((mx-mn)/101);
end

firstBin = binCtr + binSize*round( (mn-binCtr)/binSize );

tmpNbins = round((mx-binCtr)/binSize) - round((mn-binCtr)/binSize);

bins = firstBin + binSize*[0:tmpNbins];

[N, X] = hist(mtx, bins);





% RES = histoMatch(MTX, N, X)
%
% Modify elements of MTX so that normalized histogram matches that
% specified by vectors X and N, where N contains the histogram counts
% and X the histogram bin positions (see histo).

% Eero Simoncelli, 7/96.

function res = histoMatch(mtx, N, X)

if ( exist('histo') == 3 )
  [oN, oX] = histo(mtx(:), size(X(:),1));
else
  [oN, oX] = hist(mtx(:), size(X(:),1));
end

oStep = oX(2) - oX(1);
oC = [0, cumsum(oN)]/sum(oN);
oX = [oX(1)-oStep/2, oX+oStep/2];

N = N(:)';
X = X(:)';
N = N + mean(N)/(1e8);   %% HACK: no empty bins ensures nC strictly monotonic

nStep = X(2) - X(1);
nC = [0, cumsum(N)]/sum(N);
nX = [X(1)-nStep/2, X+nStep/2];

nnX = interp1(nC, nX, oC, 'linear');

if ( exist('pointOp') == 3 )
  res = pointOp(mtx, nnX, oX(1), oStep);
else
  res = reshape(interp1(oX, nnX, mtx(:)),size(mtx,1),size(mtx,2));
end


% [MIN, MAX] = range2(MTX)
%
% Compute minimum and maximum values of MTX, returning them as a 2-vector.

% Eero Simoncelli, 3/97.

function [mn, mx] = range2(mtx)

%% NOTE: THIS CODE IS NOT ACTUALLY USED! (MEX FILE IS CALLED INSTEAD)

% fprintf(1,'WARNING: You should compile the MEX version of "range2.c",\n         found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster.\n');

if (~isreal(mtx))
  error('MTX must be real-valued');  
end

mn = min(min(mtx));
mx = max(max(mtx));