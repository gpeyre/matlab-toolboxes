function [e1,e2,l1,l2] = perform_tensor_decomp(T,order)

% perform_tensor_decomp - perform an eigendecomposition.
%
%   [e1,e2,l1,l2] = perform_tensor_decomp(T);
%
%   'e1(i,j,:)' is the main eigenvector at location (i,j)
%       with associated largest eigenvalue 'l1(i,j)'.
%   'e2(i,j,:)' is the second eigenvector at location (i,j)
%       with associated smallest eigenvalue 'l2(i,j)'.
%
%   So you always have l1>=l2 (not in absolute value !).
%
%   If 'order'=='abs' then the the decomposition is done
%   so that abs(l1)>=abs(l2)
%
%   'T' must be a tensorial field (produced eg. by compute_hessian), 
%   so it should be symmetric.
%
%   Copyright (c) 2004 Gabriel Peyré

% retrieve the 4 entries of the tensor field
if ndims(T)==4  
    K11 = T(:,:,1,1);
    K12 = T(:,:,1,2);
    K21 = T(:,:,2,1);
    K22 = T(:,:,2,2);
elseif ndims(T)==2
    K11 = T(1,1);
    K12 = T(1,2);
    K21 = T(2,1);
    K22 = T(2,2);    
else
    error('T must be a tensor or a tensor field.');    
end

[n,p] = size(K11);

e1 = zeros(n,p,2);
e2 = zeros(n,p,2);
l1 = zeros(n,p);
l2 = zeros(n,p);

% trace/2
t = (K11+K22)/2;

a = K11 - t;
b = K12;

ab2 = sqrt(a.^2+b.^2);
l1 = ab2  + t;
l2 = -ab2 + t;

theta = atan2( ab2-a, b );

e1(:,:,1) = cos(theta);
e1(:,:,2) = sin(theta);
e2(:,:,1) = -sin(theta); 
e2(:,:,2) = cos(theta);

if nargin==2 & strcmp(order,'abs')
    % reorder the eigenvalue according to absolute value.
    A = abs(l1)>abs(l2);
    ee1 = prod_vf_sf(e1,A) + prod_vf_sf(e2,1-A);
    e2 = prod_vf_sf(e1,1-A) + prod_vf_sf(e2,A);
    e1 = ee1;
    ll1 = l1.*A + l2.*(1-A);
    l2 = l1.*(1-A) + l2.*A;
    l1 = ll1;
end