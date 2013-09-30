function xy = perform_classical_mds(D,ndims)


N = size(D,1);
opt.disp = 0; opt.isreal = 1; opt.issym = 1; 
M = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));
[xy, val] = eigs(M, ndims, 'LR', opt); 

for i=1:ndims
    xy(:,i) = xy(:,i)*sqrt(val(i,i));
end