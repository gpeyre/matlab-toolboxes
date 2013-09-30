function D = compute_subwindows_energy(M,k)

% compute_subwindows_matrix - compute the energy of each patch around each pixel.
%   
%   D = compute_subwindows_energy(M,k);
%
%   D(i,j) = \sum{ x \in P_{ij} }{ M(i,j)^2 }
%
%   where P_{ij} is the patch of size (2k+1)x(2k+1) around
%   pixel (i,j).
%
%   Copyright (c) 2004 Gabriel Peyré

A = symmetric_extension(M,k);
[n,p] = size(M);
D = zeros(n,p);
for x=-k:k
    for y=-k:k
        D = D + A(k+x+1:k+x+n, k+y+1:k+y+p).^2;
    end
end