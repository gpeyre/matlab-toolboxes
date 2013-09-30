function fd = div(Px,Py, options)

% div - divergence (backward difference)
%
%	fd = div(Px,Py, options);
%	fd = div(P, options);
%
%   options.bound = 'per' or 'sym'
%
%	Copyright (c) 2007 Gabriel Peyre

if size(Px,3)>1
	if nargin>1
		options = Py;
        clear Py;
	end
	Py = Px(:,:,2);
	Px = Px(:,:,1);
end

options.null = 0;
if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'sym';
end

if strcmp(bound, 'sym')
    fx = Px-Px([1 1:end-1],:);
    fy = Py-Py(:,[1 1:end-1]);
    fd = fx+fy;
else
    fx = Px-Px([end 1:end-1],:);
    fy = Py-Py(:,[end 1:end-1]);
    fd = fx+fy;
end