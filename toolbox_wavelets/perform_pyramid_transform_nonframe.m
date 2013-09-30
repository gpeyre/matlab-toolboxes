function y = perform_pyramid_transform_nonframe(x,g,options)

% perform_pyramid_transform_nonframe - fast pyramidal transform (1D or 2D).
%
%   You should use 'fwt_pyramid' instead.
%   The forward transform are the same for 'fwt_pyramid_nonframe'
%   and for 'fwt_pyramid' but the reconstruction are
%   differents. 
%   'fwt_pyramid' use the real frame operator.
%
%   For more details, see 
%   M. N. Do and M. Vetterli, Framing pyramids. Submitted IEEE Transactions on Signal Processing, December 2001
%
% y = perform_pyramid_transform_nonframe(x,g,options);
%
%   'g' is a 1D low pass filter (orthogonal or not).
%   if 'x' is vector, (forward transform), then
%       'x' is a 1D or 2D signal.
%       'y' is a cell array, y{i} is the detail at level i (1=biggest one)
%           and y{end} is the low scale decomposition.
%   if 'x' is a cell array (backward transform) then 
%       'x' is a cell array, and 'y' is a 1D or 2D signal
%
%   NB : for backward transform, the value of 
%
%   'options' is a struct that can contains
%       - 'bound': boundary extension (eiter 'sym' or 'per').
%       - 'J': number of scale for the forward transform. 
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(x)
    dir=-1;
else
    dir=1;
end
    
if dir==1
    d = ndims(x);
    if d==2 && (size(x,1)==1 || size(x,2)==1)
        d = 1;
    end
    n = length(x);
else
    d = ndims(x{1});
    if d==2 && (size(x{1},1)==1 || size(x{1},2)==1)
        d = 1;
    end
    n = length(x{1});
end

if nargin<3
    options.null = 1;
end

if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'sym';
end

if isfield(options, 'J')
    J = options.J;
else
    J = floor(log2(n));
end

% filter for second pass
if iscell(g)    % biorthogonal filters
    gg = g{2}(:);
    g = g{1}(:);    
else    % orthogonal filters
    g = g(:);
    if mod(length(g),2)==0
        gg = [g;0;0];
    else    % symmetric
        gg = g;
    end
end


if dir==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % forward transform
    
    % current low scale decomposition
    xx = x;
    
    for j=1:J
        % smooth the function
        if d==1
            % 1D convolution
            x1 = perform_convolution(xx,g, bound);
            % sub sampling
            x1 = x1(1:2:end);
            % zero padding
            x2 = zeros(size(xx));
            x2(1:2:end) = x1;
            % back convolution
            x2 = perform_convolution(x2,gg, bound);            
        else
            x1 = xx;
            % 2D : tensorial product convolution
            for i=1:size(x1,1)
                x1(i,:) = perform_convolution(x1(i,:),g, bound)';
            end
            for i=1:size(x1,2)
                x1(:,i) = perform_convolution(x1(:,i),g, bound);
            end
            % sub sampling
            x1 = x1(1:2:end, 1:2:end);
            % zero padding
            x2 = zeros(size(xx));
            x2(1:2:end, 1:2:end) = x1;
            % back convolution
            for i=1:2:size(x2,1)
                x2(i,:) = perform_convolution(x2(i,:),gg, bound)';
            end
            for i=1:size(x2,2)
                x2(:,i) = perform_convolution(x2(:,i),gg, bound);
            end
        end
        % zero padding
        y{j} = xx-x2;
        xx = x1;
    end
    
    y{end+1} = xx;
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % backward transform
    J = length(x)-1;
    
    % current low scale decomposition
    xx = x{end};
    
    for j=J:-1:1    
        if d==1
            % upsampling
            x1 = zeros(2*length(xx),1);
            x1(1:2:end) = xx;
            % back convolution
            x1 = perform_convolution(x1,gg, bound);
        else
            % upsampling
            x1 = zeros(2*size(xx));
            x1(1:2:end,1:2:end) = xx;
            % back convolution
            for i=1:2:size(x1,1)
                x1(i,:) = perform_convolution(x1(i,:),gg, bound)';
            end
            for i=1:size(x1,2)
                x1(:,i) = perform_convolution(x1(:,i),gg, bound);
            end
        end
        % adding details
        xx = x1 + x{j};
        
    end
    
    y = xx;    
    
end