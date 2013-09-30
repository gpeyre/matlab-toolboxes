function M1 = perform_image_similitude(M,u,u1,v,v1)

% perform_image_similitude
%
%   M1 = perform_image_similitude(M,u,u1,v,v1)
%
% Compute the affine similitude that map u to u1
% and v to v1, and then resample the image M.
% p and p1 are assumed to be in [0,1]²
%
%   Copyright (c) 2006 Gabriel Peyré

% T = [a -b] * [x] + [c]
%     [b  a]   [y]   [d]
%   = Q * [x;y] + t
% Solve the equations T(u)=u1 and T(v)=v1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute T %%%
% the matrix of the linear system
A = [u(1) -u(2) 1 0; ...
     u(2)  u(1) 0 1; ...
     v(1) -v(2) 1 0; ...
     v(2)  v(1) 0 1];
% the right hand size
rhs = [u1(:); v1(:)];
% solve
z = A \ rhs;
% the similitude
Q = [z(1) -z(2); z(2) z(1)];
% the translation
t = [z(3); z(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% perform resampling %%%

% original grid in the warped domain
n = size(M,1);
x = linspace(0,1,n);
[X,Y] = meshgrid(x,x);
% inverse warping P1=T^-1(P)
P = [X(:)'; Y(:)']; % position of the sampling
P(1,:) = P(1,:) - t(1); % substract translatio
P(2,:) = P(2,:) - t(2);
P1 = (Q^(-1))*P; % undo similitude
% reshape the results
X1 = reshape(P1(1,:),n,n);
Y1 = reshape(P1(2,:),n,n);

M1 = interp2(X,Y,M,X1,Y1);