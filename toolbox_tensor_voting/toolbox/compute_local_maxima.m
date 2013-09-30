function m = compute_local_maxima(M,v,d,eta)

% compute_local_maxima - find the location of local maxima.
%
%   m = compute_local_maxima(M,v,d);
%
%   M is a 2D image.
%   v is the vf indicating the direction where we should check if the 
%       function M is localy maximal (e.g. v is the gradient).
%   d is the distance used for maxima computation.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    d = 0.5;
end

v = perform_vf_normalization(v);

[n,p] = size(M);

[X,Y] = ndgrid(1:n,1:p);

Xi = X+v(:,:,1)*d;
Yi = Y+v(:,:,2)*d;
M1 = interp2(Y,X,M,Yi,Xi);

Xi = X-v(:,:,1)*d;
Yi = Y-v(:,:,2)*d;
M2 = interp2(Y,X,M,Yi,Xi);

I = find( M1<M & M2<M );
m = zeros(size(M));
m(I) = 1;