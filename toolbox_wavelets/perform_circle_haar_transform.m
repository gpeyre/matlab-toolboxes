function M1 = perform_circle_haar_transform(M, Jmin, dir, m, options)

% perform_circle_haar_transform - perform a circle-valued Haar transform.
%
%   M1 = perform_circle_haar_transform(M, Jmin, dir, m, options);
%
%   M is asumed to be known only modulo m. For instance, if m=2*pi, 
%   then M can discribes a set of angles.
%
%   If we compute 
%       M1 = perform_circle_haar_transform(M, Jmin, +1, m);
%       MM = perform_circle_haar_transform(M1, Jmin, -1, m);
%   then we have the equality MM==M modulo m, which can be tested
%   safely using
%       norm( cos(2*M-2*MM).^2 + sin(2*M-2*MM).^2, 'fro' ) is small
%
%   You can set options.normalize to make pseudo L^2-like
%   normalization of the coefficient (usefull to perform thresholding
%   after).
%   
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;
if nargin<2
    dir = 1;
end
if nargin<3
    Jmin = 1;
end
if nargin<4
    m = pi;
end

if isfield(options, 'normalize')
    normalize = options.normalize;
else
    normalize = 0;
end

if isfield(options, 'number_bands')
    number_bands = options.number_bands;
else
    number_bands = 2;
end

n = size(M,1);
Jmax = log2(n)-1;

% number of dimension
ndim = length(size(M));
if ndim==2 && ( size(M,2)==1 || size(M,1)==1 )
    ndim=1;
end

if dir==1
    M = mod(M,m);
end

if dir==-1 && normalize
    M = perform_normalization(M,Jmin,1,number_bands);
end

if size(M,1)<=2^Jmin
    M1 = M;
    return;
end

if dir==1
    M1 = fwd_step(M,m);
    if number_bands==2
        M1(1:end/2,:) = fwd_step(M1(1:end/2,:)',m)';
    else
        M1 = fwd_step(M1',m)';
    end
    options.normalize = 0;
    M1(1:end/2,1:end/2) = perform_circle_haar_transform(M1(1:end/2,1:end/2), Jmin, dir, m, options);
else
    M1 = M;
    options.normalize = 0;
    M1(1:end/2,1:end/2) = perform_circle_haar_transform(M1(1:end/2,1:end/2), Jmin, dir, m, options);
    if number_bands==2
        M1(1:end/2,:) = bwd_step(M1(1:end/2,:)',m)';
    else
        M1 = bwd_step(M1',m)';
    end
    M1 = bwd_step(M1,m);
end

if dir==1 && normalize
    M1 = perform_normalization(M1,Jmin,-1,number_bands);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M1 = fwd_step(M,m)

% coarse
C = M(1:2:end,:);
% details
D = M(2:2:end,:);

D1 = (C-D)/2;
D2 = (C+m-D)/2;
D3 = (C-(D+m))/2;

C1 = (C+D)/2;
C2 = (C+m+D)/2;
C3 = (C+D+m)/2;

% keep minimal values
v = min(abs(D1),abs(D2));
v = min(abs(v),abs(D3));
I = find(v==abs(D2));
D1(I) = D2(I);
C1(I) = C2(I);
I = find(v==abs(D3));
D1(I) = D3(I);
C1(I) = C3(I);

M1 = [mod(C1,m);D1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M1 = bwd_step(M,m)

% coarse
C = M(1:end/2,:); 
% details
D = M(end/2+1:end,:); 
M1 = M;
M1(1:2:end,:) = mod(C+D,m); 
M1(2:2:end,:) = mod(C-D,m);

function M1 = perform_normalization(M1,Jmin,s,number_bands)

ndim = 2;
n = size(M1,1);
Jmax = log2(n)-1;

% post scaling
if ndim==1
    for j=Jmin+1:Jmax
        M1(2^j:2^(j+1)) = M1(2^j:2^(j+1)) * sqrt(2^( s*(Jmax-j+1) ));
    end
    M1(1:2^Jmin) = M1(1:2^Jmin) * sqrt( 2^( s*(Jmax-Jmin+1) ) );
else
    if number_bands==2
        for j=Jmax:-1:Jmin
            a = 2^( s*(Jmax-j+1) ) / 2^(-s * 2*(Jmax-j+1));
            a = 2^( -s * 2*(Jmax-j+1));
            
            selx = 1:2^j;
            sely = (2^j+1):2^(j+1);
            M1(selx,sely) = M1(selx,sely) * a;
            
            selx = (2^j+1):2^(j+1);
            sely = 1:2^(j+1);
            M1(selx,sely) = M1(selx,sely) * (a/2);            
        end
    else
        for j=Jmax:-1:Jmin
            for q=1:3
                [selx,sely] = compute_quadrant_selection(j,q);
                M1(selx,sely) = M1(selx,sely) * 2^( s*(Jmax-j+1) );
            end
        end
    end
    M1(1:2^Jmin,1:2^Jmin) = M1(1:2^Jmin,1:2^Jmin) * 2^( s*(Jmax-Jmin+1) );
end


return;


%%% OLD CODE %%%
use_mex = 1;
if ~exist('perform_haar_transform')
    use_mex = 0;
end

if dir==-1 && use_mex
    M = reorder_coefs(M,1,-1);
end
if use_mex
    M1 = perform_haar_transform(M,Jmin, dir);
else
    M1 = perform_haar_transform_slow(M,Jmin, dir);
end
if dir==1 && use_mex
    M1 = reorder_coefs(M1,0,1);
end

if dir==1
    % post scaling
    if ndim==1
        for j=Jmin+1:Jmax
            M1(2^j:2^(j+1)) = M1(2^j:2^(j+1)) * sqrt(2^( (Jmax-j+1) ));
        end
        M1(1:2^Jmin) = M1(1:2^Jmin) * sqrt( 2^( (Jmax-Jmin+1) ) );
    else
        for j=Jmax:-1:Jmin
            for q=1:3
                [selx,sely] = compute_quadrant_selection(j,q);
                M1(selx,sely) = M1(selx,sely) * 2^( (Jmax-j+1) );
            end
        end
        M1(1:2^Jmin,1:2^Jmin) = M1(1:2^Jmin,1:2^Jmin) * 2^( (Jmax-Jmin+1) );
    end
end

% perform modulation
if 1
M1 = mod(M1,m);
tmp = min( abs(M1), abs(M1-m) );
I = find(tmp==abs(M1-m));
M1(I) = M1(I)-m;
end

if dir==1
    % post scaling
    if ndim==1
        for j=Jmin+1:Jmax
            M1(2^j:2^(j+1)) = M1(2^j:2^(j+1)) / 2^( (Jmax-j+1)/2 );
        end
        M1(1:2^Jmin) = M1(1:2^Jmin) / 2^( (Jmax-Jmin+1)/2 );
    else
        for j=Jmin:Jmax
            for q=1:3
                [selx,sely] = compute_quadrant_selection(j,q);
                M1(selx,sely) = M1(selx,sely) / 2^( (Jmax-j+1) );
            end
        end
        M1(1:2^Jmin,1:2^Jmin) = M1(1:2^Jmin,1:2^Jmin) / 2^( (Jmax-Jmin+1) );
    end
end

