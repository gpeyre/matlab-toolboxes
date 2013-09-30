function y = perform_simplex_projection_weighted(x,rho,w)
% perform_simplex_projection - compute the projection on the simplex of size rho in the diagonal metric ||.||_w
%
%   y = perform_simplex_projection_weighted(x,rho,w);
%
%   x is the projection of y on the set {a \ a >= 0 sum_i a_i = rho },
%
%   Copyright (c) 2013 Jalal Fadili
    

[n,d] = size(x);

if nargin < 3 | isempty(w)
    w = ones(n,d);
end
	
if size(w)~=[n,d]
    error('dimension of x and w must agree');
end	

if rho<=0
    error('rho should be > 0');
end

if ~isvector(rho) | numel(rho)~=n
    error('rho must be a vector of size n, x is n x d');
end

rho = repmat(rho(:),[1 d]);
	
[xs,I] = sort(x,2,'descend');
I  = I';
ws = reshape(w(sub2ind([n d],reshape(repmat([1:d],[n 1]),n*d,1),I(:))),[n d])';

xtmp = (cumsum(xs,2)-rho)./cumsum(1./ws,2);
y = max(bsxfun(@minus,x.*w,xtmp(sub2ind([n,d],(1:n)',sum(xs>xtmp./ws,2)))),0)./w;



