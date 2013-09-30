function erc = compute_erc_criterion(D,x)

% compute_erc_criterion - compute Tropp's Exact Recovery Condition
%
%   erc = compute_erc_criterion(D,x);
%
%   erc>0 ==> perfect recovery from L1 minimization.
%
%	Taken from
%		"Just relax: Convex programming methods for identifying sparse signals"
%		by J. A. Tropp.
%		IEEE Trans. Info. Theory, vol. 51, num. 3, pp. 1030-1051, Mar. 2006.
%
%   Important: the atoms of D must be unit-normed
%
%   Copyright (c) 2007 Gabriel Peyre


x = double( sign(x(:)) );
%% compute selected and un-selected columns
S1=find(abs(x)>0); % in
S2=find(x==0); % out
%% compute pseudo inverse of selected columns
D1 = D(:,S1);
warning off;
% D1 = pinv(D(:,S1));
D1 = (D1'*D1)^(-1) * D1';
% D1 = D1 ./ repmat( sqrt(sum(D1.^2,2)), [1 size(D1,2)] );
warning on;
%% compute criterion
D0 = D1 * D(:,S2);
erc = 1 - max( sum( abs(D0), 1 ) );