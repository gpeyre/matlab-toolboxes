function A = incidence2adjacency(Ic)

% incidence2adjacency - convert an incidence matrix to an adjacency matrix
%
%   A = incidence2adjacency(Ic);
%
%   If, for each edge number k of the graph linking (i,j)
%       Ic(i,k)=w1 and Ic(j,k)=w2
%   then 
%       A(i,j)=w1 and A(j,i)=w2
%
%   Copyright (c) 2006 Gabriel Peyré

%% compute list of edges
[ij,k,s] = find(sparse(Ic));
nverts = size(Ic,1);

i = ij(1:2:end);
j = ij(2:2:end);
s1 = s(1:2:end);
s2 = s(2:2:end);
i1 = [i(:); j(:)];
j1 = [j(:); i(:)];
s = [s1(:); s2(:)];

%% build sparse matrix
A = sparse(i1,j1,s,nverts,nverts);