function s = compute_skewness(x,center_mean)

% compute_skewness compute the Skewness. 
%   returns the sample skewness of the values in X.  For a
%   vector input, S is the third central moment of X, divided by the cube
%   of its standard deviation.  For a matrix input, S is a row vector
%   containing the sample skewness of each column of X.  For N-D arrays,
%   SKEWNESS operates along the first non-singleton dimension.
%
%   s = compute_skewness(x,center_mean);

% The output size for [] is a special case, handle it here.
if isequal(x,[])
    s = NaN;
    return;
end;
if nargin<2
    center_mean = 1;
end

% Figure out which dimension nanmean will work along.
sz = size(x);
dim = find(sz ~= 1, 1);
if isempty(dim)
    dim = 1;
end

% Need to tile the output of nanmean to center X.
tile = ones(1,ndims(x));
tile(dim) = sz(dim);

% Center X, compute its third and second moments, and compute the
% uncorrected skewness.
if center_mean
    x0 = x - repmat(nanmean(x), tile);
else
    x0 = x;
end
s2 = nanmean(x0.^2); % this is the biased variance estimator
m3 = nanmean(x0.^3);

denom = s2.^(1.5);
denom(abs(denom)<eps) = 1;
s = m3 ./ denom;




function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along the dimension DIM of X. 
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.13.4.2 $  $Date: 2004/01/24 09:34:32 $

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end