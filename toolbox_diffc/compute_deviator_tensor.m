function D = compute_deviator_tensor(T)

% compute_deviator_tensor - compute trace free tensor
%
%   D = compute_deviator_tensor(T);
%
%   D = T - trace(T)/2*Id
%
%   Copyright (c) Gabriel Peyre

n = size(T,1);

if size(T,3)==2 && size(T,4)==2
    t = (T(:,:,1)+T(:,:,4))/2;
    D = cat(3, t, zeros(n), zeros(n), t);
    D = T - reshape(D,[n n 2 2]);
elseif size(T,3)==3 && size(T,4)==1
    t = (T(:,:,1)+T(:,:,2))/2;
    D = cat(3, t, t, zeros(n));
    D = T - D;
else
    error('Wrong size');
end