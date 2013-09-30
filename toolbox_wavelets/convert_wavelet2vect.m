function y = convert_wavelet2vect(x, Jmin, coding_mode)

% convert_wavelet2vect - turn wavelet coefficients into a vector
%
%   y = convert_wavelet2vect(x, Jmin, coding_mode);
%
%   Turn a wavelet decomposition of an image into a vector, and reverse.
% 
%   coding_mode==1 : zig-zag ordering of each scale (respect both scale and space locality).
%   coding_mode==2 : Pack each band at a time (keep only scale locality).
%   coding_mode==3 : flush in line (destroy scale locality).
%
%   Important: for coding_mode=1/2 the scales are flushed from the finest scale
%   to coarser scales.
%
%   Copyright (c) 2005 Gabriel Peyré

n = sqrt(size(x(:),1));
Jmax = log2(n)-1;

if size(x,2)>size(x,1)
    x = x';
end

if size(x,2)>1
    if coding_mode==3
        return;
    end
    % wav2vect
    y = zeros(n^2,1); k = 1;
    for j=Jmax:-1:Jmin
        start_q = 1;
        if j==Jmin
            start_q = 0;
        end
        for q = start_q:3
            [selx,sely] = compute_quadrant_selection(j,q);
            v = x(selx,sely);
            if coding_mode==1
                % zig-zag
                [Y,X] = meshgrid(1:2^j,1:2^j);
                S = X+Y + 1e-7*(mod(X+Y,2)-0.5).*X;
                [tmp,I] = sort(S(:));
                w = v(I);
            else
                % line
                w = v(:);
            end
            y(k:k+2^(2*j)-1) = w;
            k = k + 2^(2*j);
        end
    end
else
    if coding_mode==3
        return;
    end
    % vect2wav
    for j=Jmax:-1:Jmin
        start_q = 1;
        if j==Jmin
            start_q = 0;
        end
        for q = start_q:3
            [selx,sely] = compute_quadrant_selection(j,q);
            v = x(1:2^(2*j));
            if coding_mode==1
                % zig-zag
                [Y,X] = meshgrid(1:2^j,1:2^j);
                S = X+Y + 1e-7*(mod(X+Y,2)-0.5).*X;
                [tmp,I] = sort(S(:));
                A = zeros(2^j,2^j);
                A(I) = v;
                y(selx,sely) = A;
            else
                % line
                y(selx,sely) = reshape( v ,2^j,2^j);
            end
            x(1:2^(2*j)) = [];
        end
    end
end