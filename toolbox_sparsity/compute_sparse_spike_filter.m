function h = compute_sparse_spike_filter(type,n,options)

% compute_sparse_spike_filter - compute a discrete seismic filter
%
%   h = compute_sparse_spike_filter(type,n,options);
%
%   type can be 'dergauss', 'seisfreq' or 'bump'
%   the width of the gaussian derivative is options.sigma, expressed in [-1,1]
%   so that the width in pixel is sigma*n/2.
%
%   Copyright (c) 2008 Gabriel Peyre


switch type
    case 'dergauss'
        x = linspace(-1,1,n+1)'; x(end)= [];
        sigma = getoptions(options,'sigma',.02);
        h = x.*exp( -(x.^2)/(2*sigma^2) );
        h = (1-x.^2/sigma^2).*exp( -(x.^2)/(2*sigma^2) );
        h = h-mean(h);
        %sigma = .005;
        %h = exp( -(x.^2)/(2*sigma^2) );
        h0 = h/max(h);
        h0 = h/norm(h,'fro');
        h0 = .1*h0;
        % switch for fft
        h = [h0(end/2+1:end); h0(1:end/2)];
    case {'seisfreq' 'wavelet'}
        rho = getoptions(options, 'rho', .1 * 1024/n);        
        % compute by frequency
        x = linspace(-1,1,n+1);
        hf = sin(pi*x/rho).^2;
        hf(abs(x)>rho) = 0;
        hf = [hf(n/2+1:end) hf(1:n/2-1)];
        h = real(ifft(hf));
        h = h'/max(h);
        h = h/sum(abs(h));
        
    case 'bump'
        sigma = getoptions(options,'sigma',.07);
        x = linspace(-1,1,n+1)'; x(end)=[];
        h = exp( -x.^2 / (2*(sigma)^2) );
        h = h; % /norm(h);
        h = [h(end/2+1:end); h(1:end/2)];
end
