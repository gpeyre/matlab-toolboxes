function M1 = extract_subsquares(M,w,dir)

% extract_subsquares - vectorize each sub-square of a big matrix.
%
%   M1 = extract_subsquares(M,w,dir);
%
%   M is the original matrix.
%   w is the width of the sub-squares.
%   dir=+1 for forward mapping, -1 for backward.
%
%   if M is a (p*w,q*w) matrix, then M1 is a (w^2,p*q) matrix.
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<3
    dir = 1;
end

if length(w)==1
    w = [w w];
end

if dir==1
    p = size(M);
    nb = floor( p./w );
    M1 = zeros(w(1),w(2),prod(nb));
    i = 0;
    for ky=1:nb(2)
        for kx=1:nb(1)
            i = i+1;
            selx = (kx-1)*w(1)+1:kx*w(1);
            sely = (ky-1)*w(2)+1:ky*w(2);
            M1(:,:,i) = M(selx,sely);
        end
    end
    M1 = reshape(M1, prod(w),size(M1,3));
else
    % backward transform
    % not implemented
    error('Not implemented');
end

return


%%% FAST CODE %%%
% fwd transform
n = size(M);
% sumber of sub-squares is prod(p)
p = floor(nw);
% index computation
% compute offset inside each small square
[offsy,offsx] = meshgrid(0:n:(w-1)*n,0:w-1);
offs = offsx+offsy;
offs = offs(:);
I1 = repmat(offs,1,p^2);
% compute offset due to the localisation of the square
[Y,X] = meshgrid(1:w:n,1:w:n);
I2 = X(:)' + (Y(:)'-1)*n;
I2 = repmat( I2,w^2,1 );
% add the two offsets to get the index
I = I1+I2;
% apply index reshaping
M1 = M(I);