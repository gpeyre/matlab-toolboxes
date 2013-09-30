function h = compute_histogram_rbf(f, sigma, x, options)

% compute_histogram_rbf - parzen windows density estimation
%
%   h = compute_histogram_rbf(f, sigma, x);
%
%   f is the signal, h is an estimate of the histogram, 
%   where h(i) is the density of the estimation around value x(i).
%
%   sigma is the bandwidth used to estimate the histogram (approx. size of
%       the bins).
%
%   If f is (n,2) matrix, then a joint histogram is estimated.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
    

if size(f,1)<size(f,2)
    f = f';
end
if size(x,1)<size(x,2)
    x = x';
end

nb_samples = Inf;
if isfield(options, 'nb_samples')
    nb_samples = options.nb_samples;
end
nb_samples = min(nb_samples,size(f,1));
sel = randperm(size(f,1)); sel = sel(1:nb_samples);
f = f(sel,:);

if size(f,2)==2
    h = compute_histogram_rbf_2d(f, sigma, x, options);
    return;
end

if isfield(options, 'renormalize')
    renormalize = options.renormalize;
else
    renormalize = 0;
end

n = length(x);
p = length(f);  % number of samples
f = f(:); x = x(:);

% h is of size p x n
h = repmat(f, [1 n]) - repmat( x(:)', [p 1] );
h = exp( -h.^2/(2*sigma^2) );
if renormalize
    d = repmat( sum(h,2), [1 n] );
    d(d<eps) = 1;
end
h = h./d;
h = sum(h,1);
h = h(:)/sum(h(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = compute_histogram_rbf_2d(f, sigma, x, options)

options.null = 0;

if isfield(options, 'nb_max')
    nb_max = options.nb_max;
else
    nb_max = 0;
end

if isfield(options, 'renormalize')
    renormalize = options.renormalize;
else
    renormalize = 0;
end

n = size(x,1);
p = size(f,1);

if p>nb_max
    h = zeros(n);
    niter = ceil(p/nb_max);
    for i=1:niter
        progressbar(i,niter);
        sel = (i-1)*nb_max+1:min(i*nb_max,p);
        h = h + compute_histogram_rbf_2d( f(sel,:) , sigma, x, options);
    end
    h = h/sum(h(:));
    return;
end

% hx is of size p x n x n
hx = repmat( f(:,1), [1 n n]) - repmat( reshape(x(:,1),[1 n 1]),[p 1 n] );
hy = repmat( f(:,2), [1 n n]) - repmat( reshape(x(:,2),[1 1 n]),[p n 1] );
h = exp( -(hx.^2/(2*sigma(1)^2)+hy.^2/(2*sigma(2)^2)) );
if renormalize
    d = repmat( sum(sum(h,2),3), [1 n n] );
    d(d<eps) = 1;
    h = h./d;
end
h = sum(h,1);
h = reshape(h,n,n);
h = h./sum(sum(h));