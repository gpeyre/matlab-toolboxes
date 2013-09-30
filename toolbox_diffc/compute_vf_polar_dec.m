function [theta,r] = compute_vf_polar_dec(vf)

% compute_vf_orientation - compute the polar decomposition of a vector field.
%
%   [theta,r] = compute_vf_orientation(vf)
%
%   We have vf(k1,k2,:) = [cos(theta(k1,k2));sin(theta(k1,k2))]*r(k1,k2)
%   Note that theta lies in [0,2*pi[
%
%   WORKS ONLY FOR 2D VECTOR FIELD
%
%   Copyright (c) 2004 Gabriel Peyré

if size(vf,3)~=2
    error('Works only for 2D vector fields.');
end

[theta,r] = cart2pol(vf(:,:,1),vf(:,:,2));