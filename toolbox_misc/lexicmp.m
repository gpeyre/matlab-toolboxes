function v = lexicmp(a,b, tol)

% compare using lexicographical ordering the arrays
%
%   v = lexicmp(a,b, tol);
%
%   tol (default 1e-9) is a tolerance of equality.
%
%   return -1 if a<b, +1 if a>b, 0 if a=b.
%
%   Copyright (c) 2008 Gabriel Peyre

tol = 1e-9;
a = a(:);
b = b(:);

n = max(length(a),length(b));
a(n+1:end) = -Inf;
b(n+1:end) = -Inf;

for i=1:n
    if a(i)<b(i)-tol
        v = -1; return;
    end
    if a(i)>b(i)+tol
        v = +1; return;
    end
end
v = 0;

