function y = perform_haar_transform_1d(x,dir,J, options)

% perform_haar_transform - perform the 1D haar transform along columns
%
% y = perform_haar_transform_1d(x,dir,J, options);
%
%   dir=1 for forward transform, dir=-1 for backward.
%
%   Set options.ti=1 to perform a translation invariant transform.
%
%   If x is a matrix, the transform is applied to each column.
%
%   Copyright (c) 2007 Gabriel Peyré?


s = size(x);
n = size(x,1);
m = prod(s(2:end));
x = reshape(x,[n m]);

options.null = 0;
if isfield(options, 'ti') && options.ti
    y = perform_haar_transform_ti(x,dir,J);
else
    y = perform_haar_transform_ortho(x,dir,J);
end

y = reshape(y,[size(y,1), s(2:end)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_haar_transform_ortho(x,dir,J)

n = size(x,1);
m = size(x,2);
if not(n/2^J==floor(n/2^J))
    error('n should be divisible by 2^J.');
end


if dir==1
    % x contains the coarse scale signal
    y = x;
    nj = n;
    for j=1:J
        % coarse section
        y0 = y(1:nj,:);
        y(1:nj/2,:) = (y0(1:2:end,:) + y0(2:2:end,:)) / sqrt(2);        % =y0
        y(nj/2+1:nj,:) = (y0(1:2:end,:) - y0(2:2:end,:)) / sqrt(2);     % =y1
        nj = nj/2;
    end
else
    % x contains the coarse scale signal
    nj = n/2^J;
    y = x(1:nj,:);
    for j=1:J
        % coarse section
        y1 = x(nj+1:2*nj,:); % details
        y0 = y; % coarse
        y = zeros(2*nj,m);
        y(1:2:end) = (y0+y1)/sqrt(2);
        y(2:2:end) = (y0-y1)/sqrt(2);
        nj = nj*2;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_haar_transform_ti(x,dir,J)

if dir==1
    n = size(x,1);
else
    n = size(x,1)/(J+1);
end

if dir==1
    % x contains the coarse scale signal
    y = [];
    for j=1:J
        x0 = [x(2^j+1:end,:); x(end-2^j+1:end,:)];
        % details
        y = [(x-x0)/sqrt(2); y];
        % coarse
        x = (x+x0)/sqrt(2);
    end
    % append coarse scale
    y = [x; y];
else
    % x contains the coarse scale signal
    y = x(n+1:end,:);
    x = x(1:n,:);
    for j=J:-1:1
        % retrieve details and suppress them
        d = y(1:n,:);
        y = y(n+1:end,:);
        % retrive the 2 potential coarse scales
        x0 = (x-d)/sqrt(2);
        x = (x+d)/sqrt(2);
        % average them
        x(2^j+1:end,:) = (x(2^j+1:end,:) + x0(1:end-2^j,:))/2;
    end
    y = x;
end


