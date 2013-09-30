function v = compute_crossp_vf_vf(v1,v2)

% compute_crossp_vf_vf - compute the crossproduct of 2 3D vector fields.
%
%   v = compute_crossp_vf_vf(v1,v2);
%
%   Copyright (c) 2004 Gabriel Peyre

if size(v1,3)~=3
    error('Works only for 3D vector fields.');
end

v = v1;

v(:,:,1) =  v1(:,:,2).*v2(:,:,3)-v1(:,:,3).*v2(:,:,2);
v(:,:,2) = -v1(:,:,1).*v2(:,:,3)+v1(:,:,3).*v2(:,:,1);
v(:,:,3) =  v1(:,:,1).*v2(:,:,2)-v1(:,:,2).*v2(:,:,1);