function [a,b] = compute_hodge_decompositon(v,options)

% compute_hodge_decompositon - perform Hodge decomposition of a vector field.
%
%   [a,b] = compute_hodge_decompositon(v,option);
%
%   Compute a decomposition v=a+b
%   where a is irrotational (derive from a potential)
%   and v is iccompressible (purely rotational, divergence free).
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'per';
end

d = div(v,options);
U = compute_periodic_poisson(d,strcmp(bound,'sym'));
a = grad(U,options);
b = v-a;