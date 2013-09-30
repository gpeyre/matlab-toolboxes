function A = perform_spectral_orthogonalization(A, alpha_min, alpha_max)

% perform_spectral_orthogonalization - orthogonalize a square matrix
%
% A = perform_spectral_orthogonalization(A);
%
% Copyright (c) 2006 Gabriel Peyre

if nargin<2
    alpha_min = 1;
end
if nargin<3
    alpha_max = max(1,alpha_min);
end

[U,D,V] = svd(A);

d = diag(D);
dv = abs(d);
I = find(dv>0);
du = d; du(I) = d(I) ./ abs(d(I));
dv = max(min(dv,alpha_max), alpha_min);
d = du.*dv;
D = diag(d);
% si = sign(D); D(D~=0) = 1; D = D.*si;

A = U*D*V';

return;

do_spectral = 0;
do_shuffle = 0;

if size(A,3)>1
   for i = 1:size(A,3)
       A(:,:,i) = perform_spectral_orthogonalization(A(:,:,i));
   end
   return;
end

n = size(A,1);

if do_shuffle
    A1 = A(:,2:end);
    A1 = A1(:,randperm(n-1));
    A(:,2:end) = A1;
end

if do_spectral
    [U,S,V] = svd(A);
    A = U*V';
else
    [A,R] = qr(A);
end

% enforce zero mean
A(:,1) = 1;
[A,R] = qr(A);