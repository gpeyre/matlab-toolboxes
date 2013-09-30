function D = compute_histogram_distance(H, options)

% compute_histogram_distance - compute distance between histograms
%
%   D = compute_histogram_distance(H, options);
%
%   H(:,i) is the ith histogram.
%	D(i,g) is the distance between histogram i and j.
%
%	options.histmetric is the metric used to compute the distance d(g,h) between two histograms.
%	We denote by G and H the cumulative distribution, i.e. G=cumsum(g)
%		'l2'  -> d(g,h)^2 = sum_i (g(i)-h(i))^2
%		'l1'  -> d(g,h)   = sum_i abs(g(i)-h(i))
%		'cl2' -> d(g,h)^2 = sum_i (G(i)-H(i))^2
%		'cl1' -> d(g,h)   = sum_i abs(G(i)-H(i))
%		'chi2' -> d(g,h) = sum_i (G(i)-H(i))^2 / (G(i)+H(i))
%		'bhatta' -> d(g,h) = 1 - sum_i sqrt(G(i)*H(i))
%   options.sigma is a pre-smoothing factor, counted in pixels.
%
%   Copyright (c) 2007 Gabriel Peyre

if isfield(options, 'sigma')
	sigma = options.sigma;
else
	sigma = 1.2;
end
if isfield(options, 'histmetric')
	histmetric = options.histmetric;
else
	histmetric = 'histmetric';
end

%% smooth the histograms
n = size(H,1);
m = size(H,2);
h = compute_gaussian_filter( 21,sigma/(2*n),n);
for i=1:m
    H(:,i) = perform_convolution(H(:,i),h);
end

% make them sum to 1
H = max(H,0);
H = H ./ repmat( sum(H,1), [n 1] );

%% compute distance
switch lower(histmetric)
    case {'cl2' 'cl1'}
        D = distance_matrix(cumsum(H), histmetric(2:end) );
    otherwise
        D = distance_matrix(H, histmetric );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distance_matrix(X,metric)

n = size(X,1); % dimension
p = size(X,2);

A = repmat( reshape(X', [p 1 n] ), [1 p]);
B = permute(A, [2 1 3]);

switch metric
    case 'l1'
        D = sum( abs(A-B), 3 );
    case 'l2'
        D = sqrt( sum( (A-B).^2, 3 ) );
    case 'chi2'
        a = A+B; a(a<eps) = 1;
        D = sum( ((A-B).^2) ./ a, 3 );
    case 'bhatta'
        D = 1 - sum( sqrt(A.*B), 3 );
    otherwise
        error('Unknown method');
end




