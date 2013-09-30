function V = compute_flow(p, r, s, c, n)

% compute_flow - generate a flow
% 
%   V = compute_flow(p, r, s, c, n);
%
% Generate a flow, using the equation
%           [ s  -r ]     d
%   v(x) =  [       ] * ----- + c  with d=x-p
%           [ r   s ]    |d|²
%
% If s>0 (resp. <0) then there is a source (resp. a sink)
%   at location p.
% If r>0 (resp. <0) then there is a clockwise (resp. anticlockwise)
%   vortex at location p.
%
% IMPORTANT : positions 'x' are assumed to lie in [0,1]².
%   The resulting flow is of size (n,n), regularly sampled
%   in [0,1]².
%
% If p is a (2,k) vector, then the resulting flow is 
%   the sum of the 'k' corresponding basic flow.
%   (and then 'r', 's' and 'c' should be length-k vectors).
%
%   Copyright (c) 2005 Gabriel Peyre


if nargin<5
    n= 64;
end

if size(p,1)~=2
    p = p';
end
if size(p,1)~=2
    error('p should be of size (2,k)');
end

k = size(p,2);

if nargin<4
    c = zeros(2,k);
end

if length(s)~=k || length(r)~=k
    error('r,s should be of length k.');
end
if size(c,1)~=2
    c = c';
end
if size(c,2)==1
    c = repmat(c, 1, k);
end
if size(c,1)~=2 || size(c,2)~=k
    error('c should be of size (2,k).');
end

if k>1
    V = zeros(n,n,2);
    for i=1:size(p,2)
        V = V + compute_flow(p(:,i), r(i), s(i), c(:,i), n);
    end
    return;
end

x = 0:1/(n-1):1;
[Y,X] = meshgrid(x,x);
X1 = X - p(1);
Y1 = Y - p(2);
D = X1.^2 + Y1.^2;
I = find(D<eps);
D(I) = 1;
X1 = X1./D;
Y1 = Y1./D;

V = zeros(n,n,2);
V(:,:,1) = s*X1 - r*Y1 + c(1);
V(:,:,2) = r*X1 + s*Y1 + c(2);
