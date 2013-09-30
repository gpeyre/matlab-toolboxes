function coh = compute_coherenc(D)

% compute_coherence_criterion - compute coherence criterion
%
%   [c,coh] = compute_coherence_criterion(D,x);
%
%   coh = max(|D'*D| - eye(n))
%   c = 1/2*(1+1/coh) - |x_0|
%
%   c>0 ==> perfect recovery from L1 minimization
%
%   Copyright (c) 2007 Gabriel Peyre

d = sqrt(sum(D.^2));
D = D./repmat(d, [size(D,1) 1]);

D = D'*D; D = D-diag(diag(D));
coh = max(abs(D(:)));