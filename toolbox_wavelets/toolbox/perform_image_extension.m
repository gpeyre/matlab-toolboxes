function M = perform_image_extension(M,n, type)

% perform_image_extension - extend the size by symetry
%
%   M = perform_image_extension(M,n, type);
%
%   n is the new size of the image
%
%   type can be '1side' or '2sides' depending on wether you want to extend
%   the image on only one side or on both.
%
%   Works for 2D and 3D volumetric images 
%
% Copyright (c) 2007 Gabriel Peyre

if nargin<3
    type = '1side';
end

m = size(M,1);
if strcmp(type, '1side')
    k = n-m;
    while k>size(M,1)
        M = perform_image_extension(M,size(M,1)*2,type);
        k = k - size(M,1)/2;
    end
    M = [M; M(end:-1:end-k+1,:,:)];
    M = [M, M(:,end:-1:end-k+1,:)];
    if size(M,3)>1 && size(M,3)~=3
       M = cat(3,M, M(:,:,end:-1:end-k+1));
    end
else
    w1 = floor((n-m)/2);
    w2 = n-m-w1;
    M = [M(w1:-1:1,:,:); M; M(end:-1:end-w2+1,:,:)];
    M = [M(:,w1:-1:1,:) M M(:,end:-1:end-w2+1,:)];
    if size(M,3)>1 && size(M,3)~=3
        M = cat(3, M(:,:,w1:-1:1), M, M(:,:,end:-1:end-w2+1));        
    end
end