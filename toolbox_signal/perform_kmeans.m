function [B,seeds,E] = perform_kmeans(X,nbCluster,options)

% perform_kmeans - perform the k-means clustering algorithm.
%
%   [B,seeds] = perform_kmeans(X,nbCluster,options);
%
%   'X' is a [d,n] matrix where d is the dimension of the space
%       and n is the number of points (that live in R^d).
%   'nbCluster' is the number of wanted clusters.
%
%   'B' is a vector that contain the class membership
%       of each point.
%   'seeds' is the current center of each region.
%   'E' is the energy of the segmentation (the lowest, the better).
%
%   You can provide optional parameter in a structure options: 
%     - options.seeds : You can provide initial guess for the centers of the regions
%       via a [d,nbCluster] matrix 'seeds', where seeds(:,i)
%       is the ith center point in R^d.
%     - options.intialization : if you don't provide 'seeds', then 
%       the algorithm will use a initialization depending on 'options.intialization':
%           * If options.intialization='random' : random choice of centers.
%           * If options.intialization='greedy' : greedy intialization via 
%             farthest point sampling.
%     - options.nb_iter is the number of iterations of the Lloyd algorithm.   
%
%   Copyright (c) 2004 Gabriel Peyr?


if size(X,1)>size(X,2)
    X = X';
end

if nargin<2
    nbCluster = 4;
end
options.null = 1;
if isfield(options,'nb_iter')
    nb_iter = options.nb_iter;
else
    nb_iter = 5;
end
if isfield(options,'kmeans_code')
    kmeans_code = options.kmeans_code;
else
    kmeans_code = 2;
end
if isfield(options,'etol')
    etol = options.etol;
else
    etol = 0.02;
end


if kmeans_code==1
    [B,seeds,E] = perform_kmeans_old(X,nbCluster,options);
elseif kmeans_code==2
    [B,seeds,E] = kmeansML(nbCluster,X,'maxiter',nb_iter, 'etol', etol);    
else
    [B,seeds,E] = kmeans_light(X', nbCluster, etol);  
    seeds = seeds';
end

[membership,E] = computeMembership(X,seeds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,seeds,E] = perform_kmeans_old(X,nbCluster,options);


options.null = 1;
if isfield(options,'nb_iter')
    nb_iter = options.nb_iter;
else
    nb_iter = 5;
end
if isfield(options,'initialization')
    initialization = options.initialization;
else
    initialization = 'random';
end

n = size(X,2);
d = size(X,1);  % = (2*k+1)^2

if isfield(options,'seeds') && not(isempty(options.seeds))
    seeds = options.seeds;
else
    seeds = zeros(d,nbCluster);
    if strcmp(initialization,'random')
        % original centers (selected at random)
        seed_num = floor(rand(nbCluster,1)*n)+1;
        seeds = X(:,seed_num);
    elseif strcmp(initialization,'greedy')
        % select first point at random
        seeds(:,1) = X(:,floor(rand*n)+1);
        % replace by farthest point
        D = compute_distance_to_points(X,seeds(:,1));
        [tmp,I] = max(D); I = I(1);
        seeds(:,1) = X(:,I);
        for i=2:nbCluster
            D = compute_distance_to_points(X,seeds(:,1:i-1));
            if i>2
                D = sum(D);
            end
            [tmp,I] = max(D); I = I(1);
            seeds(:,i) = X(:,I);
        end
        
    else
        error('Unknown intialization type (should be random/greedy).');
    end        
end

D = compute_distance_to_points(X,seeds);
% compute region of influence
[tmp,B] = min(D);

for i=1:nb_iter
    % compute region center
    for k=1:nbCluster
        I = find(B==k); % points belonging to cluster k
        % geometric barycenter
        Xk = X(:,I);
        if ~isempty(I)
            seeds(:,k) = sum(Xk,2)/length(I);
        else
%            warning('Empty cluster created.');
        end
    end
    % update distance to seed
    D = compute_distance_to_points(X,seeds);
    % compute region of influence
    [tmp,B] = min(D);
end

if nargout==3
   % compute the energy
   E = 0;
   for k=1:nbCluster
        I = find(B==k); % points belonging to cluster k
        D = compute_distance_to_points(X(:,I),seeds(:,k));
        E = E+sum(D);
    end
    E = E/( n*d );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary function
function D = compute_distance_to_points(X,seeds)

nbCluster = size(seeds,2);
n = size(X,2);
D = zeros(nbCluster,n);
d = size(X,1);

for k=1:nbCluster
    % distance to seed
    D(k,:) = sum( (X - repmat(seeds(:,k),1,n)).^2 );
end


% kmeans_light: K-means clustering using euclid distance.
%
%  [dataCluster codebook] = kmeans_light(data, K, stopIter)
%
%  Input and output arguments ([]'s are optional):
%   data        (matrix) of size NxD. N is the number of data (classifiee)
%                vectors,and D is the dimension of each vector.
%   K           (scalar) The number of clusters.
%   stopIter    (scalar) The threshold [0, 1] to stop learning iterations.
%                Default is .05, and smaller makes continue interations.
%   dataCluster (matrix) of size Nx1: integers indicating the cluster indicies.
%                dataCluster(i) is the cluster id for data item i.
%   codebook    (matrix) of size KxD: set of cluster centroids, VQ codewords.
%
% See also: autolabel.m
%
% Author : Naotoshi Seo
% Date   : April, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataCluster codebook distortion] = kmeans_light(data, K, stopIter)
if nargin < 3,
    stopIter = .05;
end
[N dim] = size(data);
if K > N,
    error('K must be less than or equal to the # of data');
end

% Initial codebook
codebook = data(randsample(N, K), :);

improvedRatio = Inf;
distortion = Inf;
iter = 0;
while true
    % Calculate euclidean distances between each sample and each codeword
    d = eucdist2(data, codebook);
    % Assign each sample to the nearest codeword (centroid)
    [dataNearClusterDist, dataCluster] = min(d, [], 2);
    % distortion. If centroids are unchanged, distortion is also unchanged.
    % smaller distortion is better
    old_distortion = distortion;
    distortion = mean(dataNearClusterDist);

    % If no more improved, break;
    improvedRatio = 1 - (distortion / old_distortion);
 %   fprintf('%d: improved ratio = %f\n', iter, improvedRatio); iter = iter + 1;
    if improvedRatio <= stopIter, break, end;

    % Renew codebook
    for i=1:K
        % Get the id of samples which were clusterd into cluster i.
        idx = find(dataCluster == i);
        % Calculate centroid of each cluter, and replace codebook
        codebook(i, :) = mean(data(idx, :));
    end
end

%%%% Euclidean distance matrix between row vectors in X and Y %%%%%%%
%  Input and output arguments
%   X: NxDim matrix
%   Y: PxDim matrix
%   d: distance matrix of size NxP
function d=eucdist2(X,Y);
U=~isnan(Y); Y(~U)=0;
V=~isnan(X); X(~V)=0;
d=abs(X.^2*U'+V*Y'.^2-2*X*Y');

function y = randsample(n, k, replace, w)
%RANDSAMPLE Random sample, with or without replacement.
%   Y = RANDSAMPLE(N,K) returns Y as a 1-by-K vector of values sampled
%   uniformly at random, without replacement, from the integers 1:N.
%
%   Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
%   random, without replacement, from the values in the vector POPULATION.
%
%   Y = RANDSAMPLE(...,REPLACE) returns a sample taken with replacement if
%   REPLACE is true, or without replacement if REPLACE is false (the default).
%
%   Y = RANDSAMPLE(...,true,W) returns a weighted sample, using positive
%   weights W, taken with replacement.  W is often a vector of probabilities.
%   This function does not support weighted sampling without replacement.
%
%   Example:  Generate a random sequence of the characters ACGT, with
%   replacement, according to specified probabilities.
%
%      R = randsample('ACGT',48,true,[0.15 0.35 0.35 0.15])
%
%   See also RAND, RANDPERM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $  $Date: 2003/11/01 04:28:51 $

if nargin < 2
    error('stats:randsample:TooFewInputs','Requires two input arguments.');
elseif numel(n) == 1
    population = [];
else
    population = n;
    n = numel(population);
    if length(population)~=n
       error('stats:randsample:BadPopulation','POPULATION must be a vector.');
    end
end

if nargin < 3
    replace = false;
end

if nargin < 4
    w = [];
elseif ~isempty(w)
    if length(w) ~= n
        if isempty(population)
            error('stats:randsample:InputSizeMismatch',...
                  'W must have length equal to N.');
        else
            error('stats:randsample:InputSizeMismatch',...
                  'W must have the same length as the population.');
        end
    else
        p = w(:)' / sum(w);
    end
end

switch replace
   
% Sample with replacement
case {true, 'true', 1}
    if isempty(w)
        y = ceil(n .* rand(k,1));
    else
        [dum, y] = histc(rand(k,1),[0 cumsum(p)]);
    end

% Sample without replacement
case {false, 'false', 0}
    if k > n
        if isempty(population)
            error('stats:randsample:SampleTooLarge',...
        'K must be less than or equal to N for sampling without replacement.');
        else
            error('stats:randsample:SampleTooLarge',...
                  'K must be less than or equal to the population size.');
        end
    end
    
    if isempty(w)
        % If the sample is a sizeable fraction of the population,
        % just randomize the whole population (which involves a full
        % sort of n random values), and take the first k.
        if 4*k > n
            rp = randperm(n);
            y = rp(1:k);
            
        % If the sample is a small fraction of the population, a full sort
        % is wasteful.  Repeatedly sample with replacement until there are
        % k unique values.
        else
            x = zeros(1,n); % flags
            sumx = 0;
            while sumx < k
                x(ceil(n * rand(1,k-sumx))) = 1; % sample w/replacement
                sumx = sum(x); % count how many unique elements so far
            end
            y = find(x > 0);
            y = y(randperm(k));
        end
    else
        error('stats:randsample:NoWeighting',...
              'Weighted sampling without replacement is not supported.');
    end
otherwise
    error('stats:randsample:BadReplaceValue',...
          'REPLACE must be either true or false.');
end

if ~isempty(population)
    y = population(y);
else
    y = y(:);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [membership,means,rms] = kmeansML(k,data,varargin)
% [membership,means,rms] = kmeansML(k,data,...)
%
% Multi-level kmeans.
% Tries very hard to always return k clusters.
%
% INPUT
%	k		Number of clusters
% 	data		dxn matrix of data points
%	'maxiter'	Max number of iterations. [30]
%	'dtol'		Min change in center locations. [0]
%	'etol'		Min percent change in RMS error. [0]
%	'ml'		Multi-level? [true]
%	'verbose'	Verbose level. [0]
%			    0 = none
%			    1 = textual
%			    2 = visual
%
% OUTPUT
% 	membership	1xn cluster membership vector
% 	means		dxk matrix of cluster centroids
%	rms		RMS error of model
%
% October 2002
% David R. Martin <dmartin@eecs.berkeley.edu>

% process options
maxIter = 30;
dtol = 0;
etol = 0;
ml = true;
verbose = 0;
for i = 1:2:numel(varargin),
    opt = varargin{i};
    if ~ischar(opt), error('option names not a string'); end
    if i==numel(varargin), error(sprintf('option ''%s'' has no value',opt)); end
    val = varargin{i+1};
    switch opt,
        case 'maxiter', maxIter = max(1,val);
        case 'dtol', dtol = max(0,val);
        case 'etol', etol = max(0,val);
        case 'ml', ml = val;
        case 'verbose', verbose = val;
        otherwise, error(sprintf('invalid option ''%s''',opt));
    end
end

[membership,means,rms] = ...
    kmeansInternal(k,data,maxIter,dtol,etol,ml,verbose,1);

function [membership,means,rms] = kmeansInternal(...
    k,data,maxIter,dtol,etol,ml,verbose,retry)

[d,n] = size(data);
perm = randperm(n);

% compute initial means
rate = 3;
minN = 50;
coarseN = round(n/rate);
if ~ml | coarseN < k | coarseN < minN,
    % pick random points as means
    means = data(:,perm(1:k));
else
    % recurse on random subsample to get means
    coarseData = data(:,perm(1:coarseN));
    [coarseMem,means] = ...
        kmeansInternal(k,coarseData,maxIter,dtol,etol,ml,verbose,0);
end

% Iterate.
iter = 0;
rms = inf;
if verbose>0, fwrite(2,sprintf('kmeansML: n=%d d=%d k=%d [',n,d,k)); end
while iter < maxIter,
    if verbose>0, fwrite(2,'.'); end
    iter = iter + 1;
    % Compute cluster membership and RMS error.
    rmsPrev = rms;
    [membership,rms] = computeMembership(data,means);
    % Compute new means and cluster counts.
    prevMeans = means;
    [means,counts] = computeMeans(k,data,membership);
    % The error should always decrease.
    if rms > rmsPrev, error('bug: rms > rmsPrev'); end
    % Check for convergence.
    rmsPctChange = 2 * (rmsPrev - rms) / (rmsPrev + rms + eps);
    maxMoved = sqrt(max(sum((prevMeans-means).^2)));
    if rmsPctChange <= etol & maxMoved <= dtol, 
        break; 
    end
    % Visualize.
    if verbose>1, 
        kmeansVis(data,membership,means); 
    end
end
[membership,rms] = computeMembership(data,means);
if verbose>0, fwrite(2,sprintf('] rms=%.3g\n',rms)); end

% If there's an empty cluster, then re-run kmeans.
% Retry a fixed number of times.
maxRetries = 3;
if find(counts==0),
    if retry < maxRetries,
        disp('Warning: Re-runing kmeans due to empty cluster.');
        [membership,means] = kmeansInternal( ...
            k,data,maxIter,dtol,etol,ml,verbose,retry+1);
    else
        disp('Warning: There is an empty cluster.');
    end
end

function [membership,rms] = computeMembership(data,means)
z = distSqr(data,means);
[d2,membership] = min(z,[],2);
rms = sqrt(mean(d2));

function [means,counts] = computeMeans(k,data,membership)
[d,n] = size(data);
means = zeros(d,k);
counts = zeros(1,k);
for i = 1:k,
    ind = find(membership==i);
    counts(i) = length(ind);
    means(:,i) = sum(data(:,ind),2) / max(1,counts(i));
end

%  for i = 1:n,
%    j = membership(i);
%    means(:,j) = means(:,j) + data(:,i);
%    counts(j) = counts(j) + 1;
%  end
%  for j = 1:k,
%    means(:,j) = means(:,j) / max(1,counts(j));
%  end

function z = distSqr(x,y)
% function z = distSqr(x,y)
%
% Return matrix of all-pairs squared distances between the vectors
% in the columns of x and y.
%
% INPUTS
% 	x 	dxn matrix of vectors
% 	y 	dxm matrix of vectors
%
% OUTPUTS
% 	z 	nxm matrix of squared distances
%
% This routine is faster when m<n than when m>n.
%
% David Martin <dmartin@eecs.berkeley.edu>
% March 2003

% Based on dist2.m code,
% Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

if size(x,1)~=size(y,1),
    error('size(x,1)~=size(y,1)');
end

[d,n] = size(x);
[d,m] = size(y);

% z = repmat(sum(x.^2)',1,m) ...
%     + repmat(sum(y.^2),n,1) ...
%     - 2*x'*y;

z = x'*y;
x2 = sum(x.^2)';
y2 = sum(y.^2);
for i = 1:m,
    z(:,i) = x2 + y2(i) - 2*z(:,i);
end




