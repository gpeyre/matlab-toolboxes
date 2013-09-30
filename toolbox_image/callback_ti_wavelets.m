function y = callback_ti_wavelets(x,dir,options)

% callback_atrou - callback for redundant wavelets.
%
%   y = callback_ti_wavelets(x,dir,options);
%
%   Copyright (c) 2007 Gabriel Peyre


if isfield(options, 'Jmin')
    Jmin = options.Jmin;
else
    Jmin = 4;
end

if dir==-1
    if isfield(options, 'n')
        n = options.n;
    else
        error('You must provide options.n');
    end
    Jmax = log2(n)-1;
    % convert back to a cell array
    v = {};
    for i=Jmax:-1:Jmin
        for q=1:3
            v{end+1} = reshape(x(1:n^2), n,n);
            x(1:n^2) = [];
        end
    end
    v{end+1} = reshape(x, n,n);
    x = v;
else
    n = sqrt(size(x,1));
    x = reshape(x,n,n);
end
y = perform_atrou_transform(x,Jmin,options);
if dir==1
    % convert to an array
    v = [];
    for i=1:length(y)
        v = [v; y{i}(:)];
    end
    y = v;
else
    y = y(:);
end