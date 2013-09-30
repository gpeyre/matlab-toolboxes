function [y,nbr_bits] = perform_wavelet_arithmetic_coding(x, T, options)

% perform_wavelet_arithmetic_coding - arithmetic code a wavelet transformed image.
%
% Coding:
%   [stream,nbr_bits] = perform_wavelet_arithmetic_coding(MW, T, options);
% decoding:
%   MW = perform_wavelet_arithmetic_coding(stream, T, options);
%
%	This will calls the function perform_arithmetic_coding through the scales.
%   * options.coder_type: see perform_arithmetic_coding.
%   * options.ordering_mode: specifies the ordering of coefficients
%       1 : zig-zag ordering of each scale (respect both scale and space locality).
%       2 : Pack each band at a time (keep only scale locality).
%       3 : flush in line (destroy scale locality), should not be used.
%   * options.coding_mode:
%       1 : code independantly each scale.
%       2 : code at once all the scales.
%
%   You should provide options.Jmin
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;
if isfield(options, 'coding_mode')
    coding_mode = options.coding_mode;
else
    coding_mode = 1;
end
if isfield(options, 'ordering_mode')
    ordering_mode = options.ordering_mode;
else
    ordering_mode = 1;
end
if isfield(options, 'Jmin')
    Jmin = options.Jmin;
else
    Jmin = -1;
end
if isfield(options, 'reverse_coding_order')
    reverse_coding_order = options.reverse_coding_order;
else
    reverse_coding_order = 0;
end
if nargin<2
    T = -1;
end
if size(x,1)>1 && size(x,2)>1
    dir=+1;
else
    dir=-1;
end

if dir==1 && T>0
    [tmp, x] = perform_quantization(x, T, +1);
end
if dir==1 && length(unique(x(:)))>0.5*size(x,1)^2
    warning('The data seems to be unquantized.');
    T = 1; [tmp, x] = perform_quantization(x, T);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first code the size and Jmin
nbr_bits = 0;
if dir==1
    n = size(x,1);
    % 14 bits for size
    y(1) = n; y(2) = -1;
    nbr_bits = nbr_bits + 14;
    % 3 bits for Jmin
    if Jmin<0
        warning('You should provide Jmin in options.Jmin');
        Jmin = 3;
    end
    y(3) = Jmin; y = y(:);
    nbr_bits = nbr_bits+3;
else
    n = x(1); Jmin = x(3);
    x(1:3) = [];
    y = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code coefficients
if dir==1
    x = convert_wavelet2vect(x, Jmin, ordering_mode);
end

% compute the subdivision through scales
bucket = {};
if coding_mode==1
    k = 0;
    Jmax = log2(n)-1;
    for j=Jmax:-1:Jmin
        for q=1:3
            bucket{end+1} = k+1:k+2^(2*j);
            k = k+2^(2*j);
        end
    end
    bucket{end+1} = k+1:k+2^(2*j); % low frequency
elseif coding_mode==2
    bucket{1} = 1:n^2;
else
    error('Unknown coding mode.')
end

if ~reverse_coding_order
    % code first coarser scales
    for i=1:length(bucket)
        bucket{i} = bucket{i}(end:-1:1);
    end
end

% in case of fixed arithmetic coding.
options.laplacian_type = 'genlaplacian';

% perform coding
for i=1:length(bucket)
    options.known_size = length(bucket{i});
    if dir==1
        % code
        [u,nb] = perform_arithmetic_coding(x(bucket{i}), dir, options);
        nbr_bits = nbr_bits + nb;
        % write the size 
        if i~=length(bucket)
            [y,nb] = write_integer(y,length(u(:)));
            nbr_bits = nbr_bits + nb;
        end
        % write the code
        y = [y; u(:)];
    else
        % read the size 
        if i~=length(bucket)
            [x,s] = write_integer(x);
        else
            s = length(x);
        end
        [y(bucket{i}),nb] = perform_arithmetic_coding(x(1:s), dir, options);
        x(1:s) = [];
    end
end
 
if dir==1
    %% nothing %%
else
    y = convert_wavelet2vect(y, Jmin, ordering_mode);
    if T>0
        y = perform_quantization(y, T, -1);
    end
    nbr_bits = -1;  % not defined
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,v] = write_integer(y,n)

% write_integer - write an integer into a coding stream
%
% Coding:
%   [y,nbr_bits] = write_integer(y,n);
% Decoding:
%   [y,n] = write_integer(y);

if nargin==1
    dir=-1;
elseif nargin==2
    dir=+1;
else
    error('Wrong number of arguments');
end

if dir==+1
    % code on 2 bits the range
    z = [];
    if n~=0
        while n>0
            z(end+1) = rem(n,256);
            n = fix(n/256);
        end
    else
        z = 0;
    end
    m = length(z);
    if m>4
        error('Code only 32 bits integers.');
    end
    y = [y(:); m; z(:)];
    v = 2 + 8*m; % number of bits
else
    m = y(1); y(1) = [];
    z = y(1:m); y(1:m) = [];
    v = 0;
    for i=1:m
        v = v + 256^(i-1) * z(i);
    end
end