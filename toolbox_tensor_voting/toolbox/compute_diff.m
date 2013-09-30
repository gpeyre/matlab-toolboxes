function y = compute_diff(x,order,options)

% compute_diff - compute central derivative (of order 'order') of a vector.
%
% y = compute_diff(x,order,options);
%
%   'options' is a structure:
%   - options.h is the sampling step size (default 1).
%   - options.type is the kind of finite difference.
%       type==1 is for forward differences, ie.
%           y(i) = (x(i+1)-x(i-1))/(2h), with special
%           care at boundaries.
%       type==2 is fwd differences, ie.
%           y(i) = (x(i+1)-x(i))/h, with special
%           care at boundaries.
%       type==2 is fwd differences, ie.
%           y(i) = (x(i)-x(i-1))/h, with special
%           care at boundaries.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    options.null = 0;
end
if nargin<2
    order=1;
end

if ~isfield(options, 'h')
    options.h = 1;    
end
h = options.h;

if ~isfield(options, 'type')
    options.type = 1;
end
type = options.type;


x = reshape(x, 1, prod(size(x)));

if type==1
    % central differences
    D1 = [x(2:end),x(end)];
    D2 = [x(1),x(1:end-1)];
    y = (D1-D2)/(2*h); 
    % y(1)=y(1)*2; y(end) = y(end)*2;
    y(1) = ( 4*x(2) - 3*x(1) - x(3) )/(2*h);
    y(end) = -( 4*x(end-1) - 3*x(end) - x(end-2) )/(2*h);
elseif type==2
    % fwd differences
    D1 = [x(2:end),x(end)];
    D2 = x;
    y = (D1-D2)/h;
elseif type==3
    % fwd differences
    D1 = x;
    D2 = [x(1),x(1:end-1)];
    y = (D1-D2)/h;
else
    error('This kind of differences is not supported.');
end

if order>1
    y = compute_diff(y,order-1,options);
end