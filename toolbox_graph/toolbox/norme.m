function n = norme(x)

% norme - just the L^2 norm of a vector/matrix
%
%   n = norme(x);
%
%   Copyright (c) 2003 Gabriel Peyré

if iscell(x)
    n = 0;
    for i=1:length(x)
        n = n + norme(x{i});
    end
    return;
end

n = norm(x(:),'fro');