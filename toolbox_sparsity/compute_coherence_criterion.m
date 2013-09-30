function [c,coh] = compute_coherence_criterion(D,x);

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

% normalize the atoms
coh = compute_coherence(D);
c = 1/2*(1 + 1/coh) - sum( abs(x(:))>0 );
