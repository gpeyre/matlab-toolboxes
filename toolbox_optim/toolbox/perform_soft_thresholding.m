function y = perform_soft_thresholding(x, tau)

% perform_soft_thresholding - soft thresholding
%
%   y = perform_soft_thresholding(x, tau);
%
%   y = prox_{tau*|.|_1}(x) = max(0,1-tau/|x|)*x
%
%   Proximal operator for the scalar L1 norm.
%
%   Copyright (c) 2010 Gabriel Peyre

y = max(0,1-tau./max(abs(x),1e-10)).*x;