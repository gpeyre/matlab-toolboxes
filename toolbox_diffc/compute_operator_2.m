function y = compute_operator_2(M,par,options)

% compute_operator_2 - compute a 2nd order differential
%   operator.
%
% y = compute_operator_1(M,par,options);
%
%   compute y = par(1)*d2M/dx2 + par(2)*d2M/dy2 + par(3)*d2M/dxdy 
%   where 'par' is a length-3 vector.
%
%   Copyright (c) 2004 Gabriel Peyré

H = compute_hessian(M,options);
y = par(1)*H(:,:,1,1) + par(2)*H(:,:,1,2) + par(3)*H(:,:,2,2);