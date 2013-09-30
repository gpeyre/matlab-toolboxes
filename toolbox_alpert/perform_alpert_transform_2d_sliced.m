function [w,info] = perform_alpert_transform_2d_sliced(v,pos,k,s, dir)

% perform_alpert_transform_2d_sliced - perform a forward 2D Alpert transform. You can specify a number of slice in the Y direction.
%
%   [w,info] = perform_alpert_transform_2d_sliced(v,pos,k,s, dir);
%
%   This will perform a 1.5D alpert transform on each slice  Si = { (x,y) \ (i-1)/s <= y < i/s }
%
%   'v' is the data to transform.
%   'pos' is the sampling location, a 2D vector, pos(:,i) is the ith point.
%   'k' is the number of vanishing moments (1=>Haar, 2=>linear basis ...).
%       * 'k' can be an integer, and then the algorithm will use the same
%         order 'k' for X and Y directions.
%       * 'k' can be a couple of integer k=[kx,ky] and 'kx' will be the
%         order on the X direction, and 'ky' the order on the Y direction.
%   's' is the number of slices (slices paralel to X).
%   'dir' is either 1 (fwd transform) or -1 (bwd).
%
%   'w' is the transformed data.
%   'info' is a struct containing the localisation information for each
%       basis Alpert vector.
%       'info.l' is the scale of the vector (0=coarse scale).
%       'info.n' is the space location of the vector.
%       'info.k' is the number of multiwavelet (in [1,...,k(1)*k(2)]).
%       'info.s' is the number of the slice.
%
%   Copyright (c) 2004 Gabriel Peyré

if isempty(v)
    w = [];
    return;
end

if nargin<3
    k = 2;
end
if nargin<4
    s = 1;
end
if nargin<5
    dir = 1;
end

% size of each slice, by default, equi-sized
if length(s)==1
    % special care because s can be non-integer
    s_width = ones(ceil(s),1)/s;
    ss = sum(s_width);
    s_width(end) = 1/s - (ss-1);
else
    s_width = s(:);
    s_width = s_width/sum(s_width);    % lengths should sum to 1 ...
    s = length(s);
end
s_lim = [0; cumsum(s_width)];

% number of points
P = size(pos,2);

if length(v)~=P
    error('Sampling and data does not match.');
end

w = zeros(P,1);
v = v(:);
if nargout>1
    V = zeros(P,P);
end    

x = pos(1,:);
y = pos(2,:);

if max(abs(x))<max(abs(y)/2)
    warning('It seems that your array is not properly oriented (X direction should be larger).');    
end

% scale the data
y = rescale(y);
pos(2,:) = rescale( pos(2,:) );
x = rescale(x);
pos(1,:) = rescale(pos(1,:));


info.l = ones(P,1);     % scale 
info.n = ones(P,1);     % position on X
info.k = ones(P,1);    % multiwavelet #
info.s = ones(P,1);    % number of slice


% transform each band
for i=1:ceil(s)
    y1 = s_lim(i);
    y2 = s_lim(i+1);
    % select the correct indices
    if i<s
        I = find( y>=y1 & y<y2 );
    else
        I = find( y>=y1 & y<=y2 );
    end
    % perform the transform
    [w(I),infoI] = perform_alpert_transform_2d(v(I),pos(:,I),k, dir);
    
    if isfield( infoI, 'l' )
        info.l(I) = infoI.l;     % scale 
        info.n(I) = infoI.n;     % position on X
        info.k(I) = infoI.k;     % multiwavelet #
        info.s(I) = i;           % slice #
    else
        info = [];
    end
end