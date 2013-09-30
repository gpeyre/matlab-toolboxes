function D1 = perform_floyd_algorithm(D, verbose)

% floyd - compute shortest distances on a 
%   graph using Floyd algorithm.
%
%   D1 = perform_floyd_algorithm(A, verbose);
%
%   The matrix A is the weighted adjacency matrix
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin==1
    verbose = 1;
end

N = length(D);

if verbose
    h = waitbar(0,'Computing shortest distances');
end
for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
     if verbose
         waitbar(k/N)
     end
end
if verbose
    close(h);
end

D1 = D;