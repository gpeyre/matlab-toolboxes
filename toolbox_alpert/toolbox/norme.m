function n = norme(x)

% norme - just the L^2 norm of a vector/matrix
%
%   n = norme(x);
%
%   Copyright (c) 2003 Gabriel Peyré

n = norm(x(:),'fro');