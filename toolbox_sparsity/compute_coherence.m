function coh = compute_coherence(D)

% compute_coherence - compute coherence
%
%   [coh] = compute_coherence(D);
%
%   coh = max(|D'*D| - eye(n))
%
%   Copyright (c) 2007 Gabriel Peyre

d = sqrt(sum(D.^2));
D = D./repmat(d, [size(D,1) 1]);

D = D'*D; D = D-diag(diag(D));
coh = max(abs(D(:)));