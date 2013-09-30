function werc = compute_werc_criterion(D,x,G)

% compute_werc_criterion - compute weak-erc
%
%   werc = compute_werc_criterion(D,x,G);
%
%   G is optional and is the normalized Gram matrix.
%
%   werc<1 => x can be recovered by L1 minimization.
%
%   Copyright (c) Gabriel Peyre


S = find(x~=0); % in
Sc = find(x==0); % out

% normalize dictionary
if nargin<3
    d = sqrt(sum(D.^2));
    D = D ./ repmat(d, [size(D,1),1]);
    % gram matrix with 0 on diagonal
    G = abs(D'*D);
end
G = abs(G);

G = G-diag(diag(G));
g = sum(G(:,S),2);

alpha = max(g(S));  % in
beta = max(g(Sc));  % out

werc = beta / (1-alpha);
if alpha>=1
    werc = Inf;
end

