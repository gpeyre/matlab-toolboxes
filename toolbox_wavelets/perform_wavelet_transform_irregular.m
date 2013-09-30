function y = perform_wavelet_transform_irregular(f, Jmin, dir, options)

% perform_irregular_lifting - wavelet for irregular sampling
%
%   y = perform_irregular_lifting(f, Jmin, dir, options);
%
%   Perform a wavelet transform for function sampled on points options.x.
%   For now only a linear (biorthogonal 2,2) wavelet transform is
%   implemented.
%
%   options.scaling=1 (default) rescales the coefficients
%   to mimick energy conservation.
%
%   Copyright (c) 2007 Gabriel Peyre


f = f(:);
n = size(f,1);
x = getoptions(options, 'x', linspace(0,1,n));
scaling = getoptions(options, 'scaling', 1);

Jmax = log2(n)-1;

if dir==1
    %%% FORWARD %%%

    y = [];

    for j=Jmax:-1:Jmin

        dx = x(2:end)-x(1:end-1);
        
        f1 = [f; f(end-1)];
        dx1 = [dx; dx(end)];
        % predict step
        m = ( f1(1:2:end-1).*dx1(2:2:end) + f1(3:2:end).*dx1(1:2:end-1) ) ./ ( dx1(2:2:end) + dx1(1:2:end) );
        d = f(2:2:end) - m;
        % update step
        d1 = [d(1); d];
        f = f(1:2:end) + 1/2 * ( d1(1:end-1).*dx1(2:2:end) + d1(2:end).*dx1(1:2:end-1) ) ./ ( dx1(2:2:end) + dx1(1:2:end) );
        % wavelet coefs
        if scaling
            d = d.*2^( (Jmax-j+1)/2 );
        end
        y = [d; y];
        % new sampling locations
        x = x(1:2:end);

    end
    if scaling
        f = f.*2^( (Jmax-Jmin+1)/2 );
    end
    y = [f; y];

else
    
    d0 = f; x0 = x;
    % retrieve coarse
    f = f(1:2^Jmin); d0(1:2^Jmin) = [];
    if scaling
        f = f./2^( (Jmax-Jmin+1)/2 );
    end
    for j=Jmin:Jmax
        x = x0(1:2^(Jmax-j):end);
        dx = x(2:end)-x(1:end-1);
        % retrieve details
        d = d0(1:2^j); d0(1:2^j) = [];        
        if scaling
            d = d./2^( (Jmax-j+1)/2 );
        end                
        % undo-update
        d1 = [d(1); d];
        dx1 = [dx; dx(end)];
        fe = f; f = zeros(2*size(f,1), 1);
        f(1:2:end) = fe - 1/2 * ( d1(1:end-1).*dx1(2:2:end) + d1(2:end).*dx1(1:2:end-1) ) ./ ( dx1(2:2:end) + dx1(1:2:end) );
        
        % undo predict           
        f1 = [f; f(end-1)];
        m = ( f1(1:2:end-1).*dx1(2:2:end) + f1(3:2:end).*dx1(1:2:end-1) ) ./ ( dx1(2:2:end) + dx1(1:2:end) );
        f(2:2:end) = d + m;        
    end
    y = f;

end