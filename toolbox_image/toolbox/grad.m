function [fx,fy] = grad(M, options)

% grad - gradient, forward differences
%
%   [gx,gy] = grad(M, options);
% or
%   g = grad(M, options);


options.null = 0;
if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'sym';
end

if strcmp(bound, 'sym')
    fx = M([2:end end],:)-M;
    fy = M(:,[2:end end])-M;
else
    fx = M([2:end 1],:)-M;
    fy = M(:,[2:end 1])-M;
end

if nargout==1
    fx = cat(3,fx,fy);
end