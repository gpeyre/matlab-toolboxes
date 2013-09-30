function m = mmin(A)

% mmin - minimum entry from a matrix/cell array.
%
%   m = mmin(A);
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(A)
    m = Inf;
    for i=1:length(A)
        m = min(mmin(A{i}),m);
    end
    return;
end

m = min(A(:));