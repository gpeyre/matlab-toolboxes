function x = compute_rand_sparse(n, s, model, options)

% compute_rand_sparse - compute a random sparse vector
%
%   x = compute_rand_sparse(n,s, model, options);
%
%   n is the dimension.
%   s is the number of spikes.
%   model is the distribution of the coefficients, can be 'binary', 'signed', 'gaussian', 'uniform'
%   options.nbsignals is the number of generated signal (one per column).
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    model = 'binary';
end
options.null = 0;
p = getoptions(options, 'nbsignals', 1);

x = zeros(n,p);

for i=1:p
    q = randperm(n); q = q(1:s);
    switch lower(model)
        case 'signed'
            x(q,i) = sign(randn(s,1));
        case 'binary'
            x(q,i) = 1;
        case 'gaussian'
            x(q,i) = randn(s,1);
        case 'uniform'
            x(q,i) = 2*rand(s,1)-1;
        case 'seismic'

            delta = getoptions(options, 'delta', 23); % regular spacing
            delta0 = getoptions(options, 'delta0', []); % create default that is not recovered by BP
            p = n;
            m = ceil(p/delta)+2;
            spc = [round(delta/2) ones(1,m)*delta]';
            si = sign(randn(m,1));
            k = round(length(spc)/2);
            spc(k-1:k) = delta0;
            si(k-2:k) = [1 -1 1];
            sel = cumsum(spc);
            sel(sel>p-delta/2) = [];
            si = si(1:length(sel));
            x = zeros(p,1);
            x(sel) = si;
            x = x .* (1-rand(p,1)*.3);

    end

end