function y = callback_localdct(x,dir,options)

%   callback_localdct - MCA callback for local DCT
%
% You can set options.dct_type to either:
%   'orthogonal4': overlapping orthogonal local DCT4
%   'orthogonal2': non-overlapping orthogonal local DCT2
%   'redundant': overlapping local DCT2
% The size of the windows is options.w
% The overlap is control via options.q (set to w/2 for 50% overlap)
%
%   dir=1 : forward
%   dir=-1 : backward
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
w = getoptions(options, 'w', 8);
q = getoptions(options, 'q', w/2);
n = getoptions(options, 'n', [], 1);
dct_type = getoptions(options, 'dct_type', 'redundant');
remove_lowfreq = getoptions(options, 'remove_lowfreq', 1);


if dir==+1
    tr = -1;
else
    tr = 1;
end

eta = round(w/4);
switch dct_type
    case 'redundant'
        options.normalization = 'unit';
        y = perform_windowed_dct_transform(x,w,q,n, options);
        if remove_lowfreq
            y(1:eta,1:eta,:,:) = 0;
        end
    case 'gabor'
        options.transform_type = 'fourier';
        y = perform_windowed_fourier_transform(x,w,q,n, options);
        if dir==1
            y = real(y);
        end
    case {'orthogonal4' 'orthogonal'}
        y = perform_windowed_dct4_transform(x, w, tr, options);
        if remove_lowfreq
            for dx = 0:eta-1
                for dy=0:eta-1
                    y(1+dx:w:end,1+dy:w:end) = 0;
                end
            end
        end
    case 'orthogonal2'
        y = perform_local_dct_transform(x,tr,w);
        if remove_lowfreq
            for dx = 0:eta-1
                for dy=0:eta-1
                    y(1+dx:w:end,1+dy:w:end) = 0;
                end
            end
        end
    otherwise
        error('Unknow transform');
end
