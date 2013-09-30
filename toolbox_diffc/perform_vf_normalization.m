function v = perform_vf_normalization(v)

% perform_v_normalization - renormalize a vector field.
%
%   v = perform_vf_normalization(v);
%
%   Copyright (c) 2004 Gabriel Peyre

a = nb_dims(v);
d = sqrt( sum(v.^2,a) );
d(d<1e-6) = 1;
v = v .* repmat( 1./d, [ones(a-1,1)' size(v,a)] );