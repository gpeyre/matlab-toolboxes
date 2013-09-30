function y = perform_haar_transform_slow(x, Jmin, dir);

% perform_haar_transform - compute a wavelet haar transform
%
% y = perform_haar_transform(x, Jmin, dir);
%
%   This is the non-optimized Matlab version.
%   TODO : code inverse transform (the mex version does implement it).
%
%   Copyright (c) 2004 Gabriel Peyré

% number of dimension
ndim = length(size(x));
if ndim==2 && ( size(x,2)==1 || size(x,1)==1 )
    ndim=1;
end

if ndim==1
    x = x(:);
end
n = size(x,1);
Jmax = floor( log2(n) )-1;

if nargin<3
    dir=1;
end

if dir~=1 && dir~=-1
	error('dir should be either +1 or -1.');
end


if ndim==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D transform
    if dir==1
        % x contains the coarse scale signal
        y = x;
        for j=Jmax:-1:Jmin
            n = length(x);
            n1 = ceil(n/2);
            n2 = floor(n/2);
            % fine scale
            y(n1+1:n) = ( x(1:2:2*n2) - x(2:2:2*n2) )/sqrt(2);
            % coarse scale
            if n1==n2
                y(1:n1) = ( x(1:2:n) + x(2:2:n) )/sqrt(2);
            else
                y(1:n1-1) = ( x(1:2:n-1) + x(2:2:n) )/sqrt(2);
                y(n1) = x(n);
            end
            x = y(1:n1);
        end
    else
        warning('Inverse haar transform not implemented');
        y = x;
        return;
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2D transform
    if dir==1
        for j=Jmax:-1:Jmin
            sel = 1:2^(j+1);
            x(sel,sel) = fwd_step(x(sel,sel));
            x(sel,sel) = fwd_step(x(sel,sel)')';
        end
        y = x;
    else
        for j=Jmin:Jmax
            sel = 1:2^(j+1);
            x(sel,sel) = bwd_step(x(sel,sel)')';
            x(sel,sel) = bwd_step(x(sel,sel));
        end
        y = x;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M1 = fwd_step(M)

C = M(1:2:end,:); % coarse
D = M(2:2:end,:); % details
M1 = [(C+D)/sqrt(2); (C-D)/sqrt(2)];


function M1 = bwd_step(M)

C = M(1:end/2,:); % coarse
D = M(end/2+1:end,:); % details
M1 = M;
M1(1:2:end,:) = (C+D)/sqrt(2); 
M1(2:2:end,:) = (C-D)/sqrt(2);