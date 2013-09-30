function y = reorder_coefs(x, Jmin ,dir)

% reorder_coefs - reorder the wavelet coefficients.
%
%   y = reorder_coefs(x,Jmin ,dir);
%
%   Reorder the 1D or 2D wavelet coef either
%   from 'in place' to 'ordered by scale' (if 'dir==1')
%   or from 'ordered by scale' to 'in place' (if 'dir==-1').
%
%   IMPORTANT: for the moment, works only on power of two data length.
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    error('Not enough arguments');
end
if nargin==2
    dir = +1;
end


% number of dimension
ndim = length(size(x));
if ndim==2 && ( size(x,2)==1 || size(x,1)==1 )
    ndim=1;
end

if ndim==1
    x = x(:);
end
n = size(x,1);

Jmax = log2(n)-1;

if ndim==1
    if( dir==1 )
        % from inplace to mallat
        y = [];
        for j=0:Jmax-Jmin
            % select
            sel = 1:2^j:n;  % select the coarse coef from previous step
            xs = x(sel);
            % split
            nn = length(xs);
            s = xs( 1:2:nn );  % coarse
            d = xs( 2:2:nn );  % details
            y = [d;y];
        end
        y = [s;y];
    else
        % from Mallat to inplace
        s = x;
        y = x;
        for j=0:Jmax-Jmin
            y(1:2^(j+1):end) = s(1:end/2);  % coarse
            y(1+2^j:2^(j+1):end) = s(end/2+1:end);  % coarse
            s = s(1:end/2);
        end
    end
else
    if( dir==1 )
        % from inplace to mallat
        s = x;
        y = x;
        for j=0:Jmax-Jmin
            y(1:end/2^(j+1),1:end/2^(j+1)) = s(1:2:end, 1:2:end);  % coarse/coarse
            y(1:end/2^(j+1),end/2^(j+1)+1:end/2^j) = s(1:2:end, 2:2:end);  % coarse/fine
            y(end/2^(j+1)+1:end/2^j,1:end/2^(j+1)) = s(2:2:end, 1:2:end);  % fine/coarse
            y(end/2^(j+1)+1:end/2^j,end/2^(j+1)+1:end/2^j) = s(2:2:end, 2:2:end);  % fine/fine
            s = s(1:2:end, 1:2:end);
        end
    else
        % from Mallat to inplace
        s = x;
        y = x;
        for j=0:Jmax-Jmin
            y(1:2^(j+1):end, 1:2^(j+1):end) = s(1:end/2, 1:end/2);  % coarse/coarse
            y(1:2^(j+1):end, 1+2^j:2^(j+1):end) = s(1:end/2, end/2+1:end);  % coarse/fine
            y(1+2^j:2^(j+1):end, 1:2^(j+1):end) = s(end/2+1:end, 1:end/2);  % fine/coarse
            y(1+2^j:2^(j+1):end, 1+2^j:2^(j+1):end) = s(end/2+1:end, end/2+1:end);  % fine/fine
            s = s(1:end/2, 1:end/2);
        end
    end
end