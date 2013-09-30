function y = perform_walsh_transform(x)

% perform_walsh_transform - 1D and 2D walsh transform
%
%   y = perform_walsh_transform(x);
%
%   2D transform is just tensor product of the 1D transform.
%
%   Copyright (c) 2003 Gabriel Peyré

if size(x,1)>1 && size(x,2)>1
    % 2D transform
    y = fwt2d(x);
else
    y = fwt(x);
end


function res = fwt(a)

% fwt - 1D walsh transform
%
%   y = fwt(x);
%
%   Copyright (c) 2003 Gabriel Peyré

N = length(a);
res = zeros(N,1);
if N==1
    res = a;
    return;
end

P = N/2;
a(1:P) = fwt(a(1:P));
a((P+1):N) = fwt(a((P+1):N));
for i=1:P
    tmp = a(i);
    a(i) = tmp + a(i+P);
    a(i+P) = tmp - a(i+P);
end;
res = a;

function y = fwt2d(x)

% fwt2d - 2D walsh transform
%
%   y = fwt2d(x);
%
%   Copyright (c) 2003 Gabriel Peyré

y = x; n = length(x);
for i=1:n
    y(i,:) = fwt(y(i,:)')';
end
for j=1:n
    y(:,j) = fwt(y(:,j));
end