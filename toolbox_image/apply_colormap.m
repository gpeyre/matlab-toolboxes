function A = apply_colormap(M,col)

if nargin<2
    col = 'gray';
end
switch col
    case 'gray'
        c = gray(256);
    case 'jet'
        c = jet(256);
    otherwise
        error('Unknown colormap.');
end

M = round(rescale(M,1,256));

A = zeros([size(M) 3]);
for i=1:3
    A(:,:,i) = reshape( c(M,i), size(M) );
end