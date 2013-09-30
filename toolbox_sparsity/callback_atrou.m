function y = callback_atrou(x,dir,options)

% callback_atrou - callback for redundant wavelets.
%
%   y = callback_atrou(x,dir,options);
%
%   Usefull for iterative thresholdings.
%
%   dir=1 : forward
%   dir=-1 : backward
%
%   Set options.wavtight=1 if you want the transform to be implemented as
%   tight frame.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
Jmin = getoptions(options, 'Jmin', 1);
nbdims = getoptions(options, 'nbdims', []);
tight = getoptions(options, 'wavtight', 0);

if isempty(nbdims)
    if not(iscell(x))
        nbdims = nb_dims(x);
    else
        nbdims = nb_dims(x{1});
    end
end

if dir==+1
    tr = -1;
else
    tr = 1;
end

options.ti = 1;

if nbdims==1
    y = perform_wavelet_transform(x, Jmin, tr, options);
else
    if iscell(x) && tight
        %% normalize
        for i=1:(length(x)-1)/3
            x{3*i-2} = x{3*i-2}*2^(i);
            x{3*i-1} = x{3*i-1}*2^(i);
            x{3*i  } = x{3*i  }*2^(i);
        end
        x{end} = x{end}*2^(i);
    end
    y = perform_wavelet_transform(x,Jmin,tr, options);
    if not(iscell(x)) && tight
        %% normalize
        for i=1:(length(y)-1)/3
            y{3*i-2} = y{3*i-2}/2^(i);
            y{3*i-1} = y{3*i-1}/2^(i);
            y{3*i  } = y{3*i  }/2^(i);
        end
        y{end} = y{end}/2^(i);
    end
end

if nbdims==1 && dir==-1
    y = y(:);
end


function d = nb_dims(x)

% nb_dims - debugged version of ndims.
%
%   d = nb_dims(x);
%
%   Copyright (c) 2004 Gabriel Peyré

if isempty(x)
    d = 0;
    return;
end

d = ndims(x);

if d==2 && (size(x,1)==1 || size(x,2)==1)
    d = 1;
end