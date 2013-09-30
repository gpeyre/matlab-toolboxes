function M = perform_local_dct_transform(M,dir,w)

% perform_local_dct_transform - perform a JPEG-like transfrom
%
%   M = perform_local_dct_transform(M,dir,w, q);
%
%   dir=+1 : forward transform
%   dir=-1 : backward transform
%
%   w is the size of the blocks.
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<3
    w = 8;
end

n = size(M,1);


for i=1:floor(n/w);
    for j=1:floor(n/w);
        selx = (i-1)*w+1:min(i*w,n);
        sely = (j-1)*w+1:min(j*w,n);
        if dir==1
            M(selx,sely) = perform_dct_transform(M(selx,sely),+1);
        else
            M(selx,sely) = perform_dct_transform(M(selx,sely),-1);
        end
    end
end
