
n = 100;  %  number of measurements
m = n*10;
s = 10;

D = compute_compressed_sensing_matrix(n,m);

% random signal
x = compute_rand_sparse(m,s);
% sensed signal
y = D*x;

niter_extension = 30; 
fuchs_earlystop = 0;

% original fuchs
[f,i,si] = compute_fuchs_criterion(D,x);

return; 

% extended fuchs
[x1,succeded,flist] = perform_extension_support(D,x,niter_extension,fuchs_earlystop);
fe = min(flist);

x1 = perform_l1_recovery(D,y,options);