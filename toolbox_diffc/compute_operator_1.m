function y = compute_operator_1(M,par,options)

% compute_operator_1 - compute a 1st order differential
%   operator.
%
% y = compute_operator_1(M,par,options);
%
%   compute y = <grad(M),par> where 'par' is a length-2 vector.
%
%   Copyright (c) 2004 Gabriel Peyré

grad = compute_grad(M,options);
y = par(1)*grad(:,:,1)+par(2)*grad(:,:,2);