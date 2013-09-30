function T = compute_tensor_field(n, pos,type, options)

% compute_tensor_field - compute a parametric 2D tensor field
%
%   T = compute_tensor_field(n, pos,type, options);
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if iscell(pos)
    T = zeros(n,n,3);
    for i=1:min(length(pos),length(type))
        T = T + compute_tensor_field(n, pos{i},type{i}, options);
    end
    return;
end

s = getoptions(options, 'tensor_sigma', .1);
bound = getoptions(options, 'bound', 'sym');

x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);
X = X-pos(1);
Y = Y-pos(2);

if strcmp(bound, 'per')
    X(X>1/2) = X(X>1/2) - 1; X(X<-1/2) = X(X<-1/2) + 1;
    Y(Y>1/2) = Y(Y>1/2) - 1; Y(Y<-1/2) = Y(Y<-1/2) + 1;
end
d = exp( -(X.^2+Y.^2) / (2*s^2) );

switch type
    case 'wedge'
        A = cat(3, X,-X,Y);
    case 'trissector'
        A = cat(3, X,-X,-Y);        
    case 'node'
        A = cat(3, X.^2-Y.^2,X.^2-Y.^2,2*X.*Y);
    case 'center'
        A = cat(3, Y.^2-X.^2,X.^2-Y.^2,-2*X.*Y);
    case 'saddle'
        A = cat(3, X.^2-Y.^2,X.^2-Y.^2,-2*X.*Y);
    otherwise
        error('Unknown singularity type.');
end
T = cat(3,d,d,d) .* A;