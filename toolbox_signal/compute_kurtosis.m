function k = compute_kurtosis(x, center_mean)
%compute_kurtosis - compute the Kurtosis. 
%
%   returns the sample kurtosis of the values in X.  For a
%   vector input, K is the fourth central moment of X, divided by fourth
%   power of its standard deviation.  For a matrix input, K is a row vector
%   containing the sample kurtosis of each column of X.  For N-D arrays,
%   KURTOSIS operates along the first non-singleton dimension.
%
%   K = compute_kurtosis(X, center_mean);
%
%   Set center_mean=0 if you do not want to remove the mean before
%   computing the kurtosis.

% The output size for [] is a special case, handle it here.
if isequal(x,[])
    k = NaN;
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

% Center X, compute its fourth and second moments, and compute the
% uncorrected kurtosis.
if center_mean
    x0 = x - repmat(nanmean(x), tile);
else
    x0 = x;
end

s2 = nanmean(x0.^2); % this is the biased variance estimator
m4 = nanmean(x0.^4);


denom = s2.^2;
denom(abs(denom)<eps) = 1;
k = m4 ./ denom;


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

