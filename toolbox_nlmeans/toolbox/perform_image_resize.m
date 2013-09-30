function newimg = image_resize(img,p1,q1)

% im_resize - resize an image using bicubic interpolation
%
%   newimg = im_resize(img,nw,nh);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin ~= 3
    error('usage: im_resize(image,new_wid,new_ht)');
end;

p = size(img,1);
q = size(img,2);

[X,Y] = meshgrid( (0:p-1)/(p-1), (0:q-1)/(q-1) );
[XI,YI] = meshgrid( (0:p1-1)/(p1-1) , (0:q1-1)/(q1-1) );


newimg = interp2( X,Y, img, XI,YI ,'cubic');