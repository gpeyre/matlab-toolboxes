function X = compute_bezier_curve(C,N)

% compute_bezier_curve - plots N points on the Bezier curve with control points in C
%
%   compute_bezier_curve(C,N);
%
%   C is a matrix of size (2,n)
%   X is a matric of size (2,N)
%
%   The bezier curve is of degree n-2.
%   If no output, then the curve is displayed together with the control
%   points.
%
%   Copyright (c) 2006 Gabriel Peyré

if size(C,2)<size(C,1)
    C = C';
end
if size(C,1)~=2
    error('Works only for 2D curves');
end

n = length(C(1,:));            % number of control points
a = linspace(0,1,N);             % parameter values along the Bezier curve
I = (0:(n-1))'*ones(1,N);      % matrix of index values
P = ones(n,1)*a;               % matrix of probabilities
A = binopdf(I,n-1,P);          % weights for control points
X = C*A;                       % points on the Bezier curve

if nargout==0
    hold on;
    plot(X(1,:), X(2,:));
    plot(C(1,:), C(2,:), '.');
    hold off;
end