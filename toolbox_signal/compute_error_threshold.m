function [ffM,M] = compute_error_threshold(y,T);

% compute_error_threshold - compute the error for an approximation in an orthogonal basis.
%
%   [ffM,M] = compute_error_threshold(y,T);
%
%   'y' are the coefficient of the approximation in the base.
%   'T' is the threshold for the approximation.
%
%   Copyright (c) 2004 Gabriel Peyré

I = find(abs(y)>=T);
y(I) = 0;
ffM = norme(y);
M = length(I);