function M2 = perform_zooming(M,q)

n = size(M,1);
s = size(M,3);

M1 = M;
M1 = crop(M1,q);
x = linspace(0,1,q);
xi = linspace(0,1,n);
[X,Y] = meshgrid(x,x);
[Xi,Yi] = meshgrid(xi,xi);
M2 = M;
for j=1:s
    M2(:,:,j) = interp2(X,Y,M1(:,:,j),Xi,Yi);
end