function T = perform_tensor_mapping(T,dir)

% perform_tensor_mapping - go back and forth from tensor to non-linear domain
%
% forward:
%   U = perform_tensor_mapping(T,+1);
% backward:
%   T = perform_tensor_mapping(U,-1);
%
%   T is an (n,n,3) tensor field, with egenvalues L1 and L2.
%   U(:,:,1) is the energy (L1+L2)
%   U(:,:,2) is the anisotropy ((L1-L2)/(L1+L2))^2
%   U(:,:,3) is the orientation (angle of the main eigenvector.
%
%   Works only for 2D fields.
%
%   See also: perform_tensor_decomp.
%
%   Copyright (c) 2007 Gabriel Peyre

if iscell(T)
    for i=1:length(T)
        T{i} = perform_tensor_mapping(T{i},dir);
    end
    return;
end


if dir==1
    if size(T,3)==3 && size(T,4)==1
        a = T(:,:,1);
        b = T(:,:,2);
        c = T(:,:,3);
    elseif size(T,3)==2 && size(T,4)==2
        a = T(:,:,1,1); 
        b = T(:,:,2,2);
        c = T(:,:,1,2);
    else
        error('Tensor field has wrong size');        
    end
    % energy
    T(:,:,1) = a+b;
    % anisotropy
    d = (a+b).^2; d = max(d, 1e-9);
    T(:,:,2) = 4*(a.*b-c.^2)./d;
    % orientation
    T(:,:,3) = 1/2 * atan2( 2*c, (a-b) );
else
    E = T(:,:,1);
    A = T(:,:,2);
    O = T(:,:,3);
    P = 1/4 * A .* E.^2; %  product (determinant)
    D = E.^2-4.*P; D(D<0) = 0;
    l1 = (E+sqrt(D)) / 2;
    l2 = E-l1;
    e1 = cat(3, cos(O), sin(O));
    e2 = cat(3, -sin(O), cos(O));
    T = perform_tensor_eigendecomposition(e1,e2,l1,l2);
end