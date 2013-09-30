function H = perform_square_extraction(M,epsilon,w)

% forward
%   H = perform_square_extraction(M,epsilon,w);
% backward
%   M = perform_square_extraction(H,epsilon,n);
%
%   Copyright (c) 2006 Gabriel Peyré

if size(M,4)==1
    dir = 1;
    n = size(M,1);
else
    dir = -1;
    n = w;
    w = size(M,1);
end

n1 = ceil(n/w)*w;

% sampling locations
[dY,dX] = meshgrid(0:w-1,0:w-1);
[Y,X] = meshgrid( (1:w:n1)+epsilon(2), (1:w:n1)+epsilon(1) );
m = length(X(:));
X = reshape(X(:),[1 1 1 m]);
Y = reshape(Y(:),[1 1 1 m]);
IX = repmat(dX,[1 1 1 m]) + repmat( X, [w w 1 1] );
IY = repmat(dY,[1 1 1 m]) + repmat( Y, [w w 1 1] );
IX = mod(IX-1,n)+1;
IY = mod(IY-1,n)+1;
I = sub2ind([n n], IX, IY);


for s=1:size(M,3)
    if dir==1
        U = M(:,:,s);
        H(:,:,s,:) = U(I);
    else
        U(I) = M(:,:,s,:);
        H(:,:,s) = reshape(U,n,n);
    end
end