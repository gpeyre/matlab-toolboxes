function m = mmax(A)

% mmax - maximum entry from a matrix.
%
%   m = mmax(A);
%
%   Works for matrix and cell array.
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(A)
    m = -Inf;
    for i=1:length(A)
        m = max(mmax(A{i}),m);
    end
    return;
end

m = max(A(:));