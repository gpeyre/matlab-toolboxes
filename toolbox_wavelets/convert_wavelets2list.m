function y = convert_wavelets2list(x, Jmin)

% convert_wavelets2list - convert a wavelet transform into a cell array.
%
%	y = convert_wavelets2list(x, Jmin);
%
%   If x is a 2D array, then y is a cell array.
%       y{3*(j-1)+q} contains wavelets coefficients at scale j and orientation q.
%   If x is a cell array, then y is an image.
%
%   Copyright (c) 2007 Gabriel Peyre


if not(iscell(x))
    y = {};
    n = size(x,1);
    Jmax = log2(n)-1;
    
else
    Jmax = log2(size(x{1},1));
    n = 2^(Jmax+1);
    y = zeros(n);    
end

m = 0;
for j=Jmax:-1:Jmin
   qmin = 1-(j==Jmin);
   for q = qmin:3
        m = m+1;
        [selx,sely] = compute_quadrant_selection(j,q);
        if not(iscell(x))
            y{end+1} = x(selx,sely);
        else
            y(selx,sely) = x{m};
        end
       
   end
end
    



return;


if nargin<4
    pack_imagettes = 1;
end

if ~iscell(x)
    d = size(x,3);
    if nargin<3
        Jmax = log2(size(x,1))-1;
    end
    if nargin<2
        error('You must provide Jmin');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % turn the wavelets coefficient into a cell array
    y = {};
    for j=Jmax:-1:Jmin
        % concatenate the 3 imagettes
        M_list = {};
        for q=1:3
            for m=1:d
                [selx,sely] = compute_quadrant_selection(j,q);
                if pack_imagettes
                    M_list{end+1} = x(selx,sely,m);
                else
                    y{end+1} = x(selx,sely,m);
                end
            end
            if ~pack_imagettes
                y = {y{:} M_list{:}};
                M_list = {};
            end
        end
        if pack_imagettes
            y{end+1} = M_list;
        end
    end
    % add the low scale image
    M_list = {};
    for m=1:d
        if pack_imagettes
            M_list{end+1} = x(1:2^Jmin,1:2^Jmin,m);
        else
            y{end+1} = x(1:2^Jmin,1:2^Jmin,m);
        end
    end
    if pack_imagettes
        y{end+1} = M_list;
    end
else
    if pack_imagettes
        Jmax = log2(size(x{1}{1},1));
        d = length(x{1})/3;
        Jmin = Jmax - length(x)+2;
    else
        Jmax = log2(size(x{1},1));
        d = length(x{1});
        Jmin = Jmax - (length(x)-1)/3+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % turn the cell array into a list
    n = 2^(Jmax+1);
    y = zeros(n,n,d);
    k=0;
    for j=Jmax:-1:Jmin
        if pack_imagettes
            k = k+1;
            M_list = x{k};
        end
        % concatenate the 3 imagettes
        for q=1:3
            if ~pack_imagettes
                k = k+1;
                M_list = {x{k}};
            end
            for m=1:d
                [selx,sely] = compute_quadrant_selection(j,q);
                if pack_imagettes
                    y(selx,sely,m) = M_list{q+(m-1)*d};
                else
                    y(selx,sely,m) = M_list{m};
                end
            end
        end
    end
    % add the low scale image
    if pack_imagettes
        M_list = x{end};
    else
        M_list = x{end-d+1:end}
    end
    for m=1:d
        y(1:2^Jmin,1:2^Jmin,m) = M_list{m};
    end
end