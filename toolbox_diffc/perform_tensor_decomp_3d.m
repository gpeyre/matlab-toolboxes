function [out1,out2,out3] = perform_tensor_decomp_3d(in1,in2)

% perform_tensor_decomp_3d - decompose a 3D tensor field
%
% Decomposition:
%   [U,Lambda,err] = perform_tensor_decomp_3d(T);
% Re-composition:
%   T = perform_tensor_decomp_3d(U,Lambda);
%
%   Copyright (c) Gabriel Peyre 2008

if nargin==2
    % re-composition
    U = in1; Lambda = in2;
    d = size(U);
    U = reshape(U, [prod(d(1:end-2)) d(end-1) d(end)]);
    Lambda = reshape(Lambda, [prod(d(1:end-2)) d(end-1)]);    
    L = U*0; 
    L(:,1,1) = Lambda(:,1);
    L(:,2,2) = Lambda(:,2);
    L(:,3,3) = Lambda(:,3);
	out1 = prod_tf_tf( U, prod_tf_tf( L, permute(U, [1 3 2]) ) );
    out1 = reshape(out1, d);
    return;
end

T = in1;
d = size(T);
T = reshape(T, [prod(d(1:end-2)) d(end-1) d(end)]);
n = size(T,1);
if d(end-1)~=3 || d(end)~=3
    error('Works only for 3D fields');
end

% U(:,:,i) is the ith eigenvector
U = randn(n,3,3);

% find first eigenvector field
niter = 100;
v = randn(n,3);
for iter=1:niter
    % v <- T*v/|T*v|
    v = perform_vf_normalization( prod_tf_vf(T,v) );
end
U(:, :,1) = v;

% orthogonal basis of the orthogonal of the eigenvector
% use cross product
U(:,:,3) = cross( U(:,:,1),U(:,:,2), 2 );
U(:,:,3) = perform_vf_normalization( U(:,:,3) );
U(:,:,2) = cross( U(:,:,1),U(:,:,3), 2 ); 
U(:,:,2) = perform_vf_normalization( U(:,:,2) );

% change of basis T1=U'*T*U
T1 = prod_tf_tf( permute(U,[1 3 2]), prod_tf_tf( T, U ) );

% apply diagonalization to sub-tensor
TT = reshape( T1(:,2:3,2:3), [n 1 2 2] );
[e1,e2] = perform_tensor_decomp( TT );
e1 = squeeze(e1); e2 = squeeze(e2);
U1 = zeros(n, 3,3);
U1(:,1,1) = 1;
U1(:,2:3,2) = e1;
U1(:,2:3,3) = e2;

% compute U with U1 to get the final orthogonal basis
U2 = prod_tf_tf(U,U1);

% perform diagonalization to retrieve eigenvalue
L = prod_tf_tf( permute(U2,[1 3 2]), prod_tf_tf( T, U2 ) );
Lambda = cat(4, L(:,1,1), L(:,2,2), L(:,3,3) );
err = sum( abs(L(:,1,2))+abs(L(:,1,3))+abs(L(:,2,3)) ) / (3*n);

U2 = reshape(U2, d);
Lambda = reshape(Lambda, d(1:end-1));

out1 = U2;
out2 = Lambda;
out3 = err;