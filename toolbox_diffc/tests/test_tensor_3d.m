% test for 3D tensor diagonalization

% boundary conditions
options.bound  = 'sym';
% order of the derivative
options.order = 2;  % centered differences

n = 30;
M = randn(n,n,n);
% compute hessian tensors
T = compute_hessian(M,options);


[U,Lambda,err] = perform_tensor_decomp_3d(T);
T1 = perform_tensor_decomp_3d(U,Lambda);

% reconstruction error, should be small
sqrt( mean( (T(:)-T1(:)).^2 ) )
