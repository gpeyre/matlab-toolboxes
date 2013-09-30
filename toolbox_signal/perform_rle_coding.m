function [y1,y2,histos] = perform_rle_coding(x, dir, options)

% perform_rle_coding - perform run length coding of a 1D sequence
%
% Coding: 
%   [stream,nb_bits] = perform_rle_coding(x, +1, options);
% De-Coding: 
%   x = perform_rle_coding(stream, +1, options);
%
%   Warning : works only for binary 0/1 data.
%
%   options.rle_coding_mode == 'shannon' : only provides an estimate of nb_bits
%   options.rle_coding_mode == 'arithmetic' : perform coding using artithmetic coder
%   options.rle_coding_mode == 'arithfixed' : uses laplacian density
%       estimation for coding (usually works better).
%
%   Copyright (c) 2006 Gabriel Peyré


options.null = 0;
if isfield(options, 'rle_coding_mode')
    rle_coding_mode = options.rle_coding_mode;
else
    rle_coding_mode = 'shannon';
end

if isempty(x)
    y1 = [];
    y2 = 0;
    histos = [];
    return;
end
    
if dir==1
    
    if size(x,1)==1 && size(x,2)>1
        x = x(:);
    end
    
    % force to binary
    x = rescale(x);
    x = x>0.5;
    n = size(x,1); % length of signals
    
    if size(x,2)>1
        error('Do not work for multiple signals.');
    end

    % compute the run-lengths    
    i = [ find(x(1:end-1) ~= x(2:end));  n ];
    r = diff([ 0; i ]);
    rv = x(i);
    R{1} = r(1:2:end);
    R{2} = r(2:2:end);
            
    % code the runs of 1st symbol and then the runs of second symbol
    nb_bits = 1;
    streams = {[],[]};

    for i=1:2
        r = R{i};
        stream = [];
        switch rle_coding_mode
            
            case 'arithmetic'
                %% arithemtic coding of the sequence %%
                opt.coder_type = 1;
                if ~isempty(r)
                    [stream,nb] = perform_arithmetic_coding(r, +1, opt);
                else
                    nb = 1;
                end
                nb = nb - 12; % we know the bound and the size
                if isfield(options, 'w')
                    nb3 = log2(prod(options.w)) * length(r);
                    stream = r;
                    if nb3<nb
                        nb = nb3;
                    end
                end
                
            case 'arithfixed'
                nb = 0;
                % code with generalized lapalcian fixed encoding
                opt.coder_type = 7; % generalized laplacian coding
                opt.laplacian_type = 'genlaplacian0'; % positive values
                [st,nb] = perform_arithmetic_coding(r, +1, opt);
                stream(1) = 1;
                % to compare with classical arithmetic encoding
                [st2,nb2] = perform_arithmetic_coding(r, +1);
                if (nb2<nb) & 0
                    nb = nb2;
                    st = st2;
                    stream(1) = 2;
                end
                % to compare with brute force encoding
                if 0  % isfield(options, 'w') && 0
                    nb3 = log2(prod(options.w)) * length(r);
                    st = r;
                    if nb3<nb
                        nb = nb3;
                        stream(1) = 3;
                    end
                end
                % add the code to choose between the coders
                nb = nb+1;
                stream = [stream; st];
                
            case 'shannon'
                m = max(r);
                %% shannon estimation %%
                hn = compute_histogram(r) * length(r);
                % then retrieve the histograms used for coding
                if isfield(options, 'histos')
                    h = options.histos{i};
                    h(length(h)+1:m) = 1e-16; h = h/sum(h);
                else
                    % use self histograms (works if the signal is really long)
                    if length(r)>0
                        h = hn/length(r);
                    else
                        h = hn;
                    end
                end
                I = find(h==0); h(I) = eps;
                histos{i} = h;
                nb = - sum( hn .* log2(h(1:m)) );
                
            case 'nocoding'
                % send directly the integers
                stream = [stream(:); length(r); r(:)];
                % estimate coding cost
                nb = log2(max(r)) * length(r);
                
            otherwise
                error('Unknown coding');
        end
        nb_bits = nb_bits + nb;
        streams{i} = stream;
    end
    
    % multiplex the two signals and add the intial bit
    stream = [x(1)];
    p = length(streams{1});
    [stream,nb] = append_integer(stream, p);
    stream = [stream; streams{1}; streams{2}]; 
    nb_bits = nb_bits + nb;
    % output
    y1 = stream; y2 = nb_bits;
    
else
    
    %%% decoding %%%
    stream = x; R = {};
    
    % retrieve first bit
    x0 = stream(1);  stream(1) = [];
    % retrieve length of first signal
    [stream,p] = append_integer(stream);
    % retrieve the two flows of bits
    streams = {stream(1:p), stream(p+1:end)};
    
    for i=1:2
        stream = streams{i};
        switch rle_coding_mode
            case 'arithmetic'
                
                error('Not implemented.');
                
            case 'arithfixed'
                
                type_code = stream(1); stream(1) = [];
                
                switch type_code
                    case 1                        
                        opt.coder_type = 7; % generalized laplacian coding
                        opt.laplacian_type = 'genlaplacian0'; % positive values
                        r = perform_arithmetic_coding(stream, -1, opt);
                    
                    case 2
                        
                        % to compare with classical arithmetic encoding
                        r = perform_arithmetic_coding(stream, -1);
                    
                    otherwise 
                        error('Problem with decoding.');
                end
                
            case 'nocoding'
                s = stream(1); stream(1) = [];
                r = stream(1:s); stream(1:s) = [];
                
            otherwise
                error('Unknown coding');
        end
        R{i} = r;
    end
    x = [];
    %%% recompose everything %%%
    k = 0;
    while ~isempty(R{k+1})
        x = [x; ones(R{k+1}(1),1)*k];
        R{k+1}(1) = [];
        % switch of bit
        k = 1-k;
    end
    if x0==1
        x = 1-x;
    end
    if ~isempty(R{1}) || ~isempty(R{2})
        error('Problem with decoding.');
    end
    y1 = x; y2 = -1;
end


function [stream,v] = append_integer(stream,n)

% append_integer - append an integer at the end of a coding stream
%
% Coding:
%   [stream,nbr_bits] = append_integer(stream,n);
% Decoding:
%   [stream,n] = append_integer(stream);
%
%   The integer is assumed to be positive.
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin==1
    dir=-1;
elseif nargin==2
    dir=+1;
else
    error('Wrong number of arguments');
end

if dir==+1
    if n<0
        error('Code only positive integer.');
    end
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
    stream = [stream(:); m; z(:)];
    v = 2 + 8*m; % number of bits
else
    m = stream(1); stream(1) = [];
    z = stream(1:m); stream(1:m) = [];
    v = 0;
    for i=1:m
        v = v + 256^(i-1) * z(i);
    end
end

function [stream,y] = append_stream(stream, st)

% append_stream - add a stream at the end of another one
%
% Coding:
%   [stream,nb_bits] = append_stream(stream, st)
% Decoding
%   [stream,st] = append_stream(stream)
%
%   nb_bits is the additional coding cost.
%
% Copyright (c) 2006 Gabriel Peyré



if nargin==1
    dir=-1;
elseif nargin==2
    dir=+1;
else
    error('Wrong number of arguments');
end

if dir==+1
    % write size
    [stream,nb_bits] = append_integer(stream, length(st));
    % write stream
    stream = [stream(:); st(:)];
    y = nb_bits;
else
    % retrieve size
    [stream,nb] = append_integer(stream);
    % read stream
    st = stream(1:nb); stream(1:nb) = [];
    y = st;
end