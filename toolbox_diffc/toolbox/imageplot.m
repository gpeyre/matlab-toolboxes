function h1 = imageplot(M,str, a,b,c)

% imageplot - diplay an image and a title
%
% Example of usages:
%   imageplot(M);
%   imageplot(M,title);
%   imageplot(M,title,1,2,1);   % to make subplot(1,2,1);
%
%   If you want to display several images:
%       imageplot({M1 M2}, {'title1', 'title2'});
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<2
    str = [];
end

nbdims = nb_dims(M);

if iscell(M)
    q = length(M);
    if nargin<5
        c = 1;
        a = ceil(q/4);
        b = ceil(q/a);
    end
    if (c-1+q)>(a*b)
        warning('a and c parameters not large enough');
        a = ceil((c-1+q)/4);
        b = ceil((c-1+q)/a);
    end
    for i=1:q
        if iscell(str)
            str1 = str{i};
        else
            str1 = str;
        end
        h{i} = imageplot(M{i},str1, a,b,c-1+i);
    end
    global axlist;
    if not(isempty(axlist))
        linkaxes(axlist, 'xy');
    end

    if nargout>0
        if exist('h')
            h1 = h;
        else
            h1 = [];
        end
    end
    return;
end

if nargin==5
    global axlist;
    global imageplot_size;
    if c==1 || isempty(imageplot_size) || imageplot_size~=size(M,1)
        clear axlist; 
        global axlist; 
        axlist = [];
        imageplot_size = size(M,1);
    end
    axlist(end+1) = subplot(a,b,c);
end

if nbdims==1
    h = plot(M); axis tight;
elseif size(M,3)<=3
    % gray-scale or color image
    if size(M,3)==2
        M = cat(3,M, zeros(size(M,1),size(M,2)));
    end
    if size(M,3)==1
        colormap gray(256);
    else
        colormap jet(256);
    end
    h = imagesc(rescale(M)); axis image; axis off;
else
    % volumetric image
    h = vol3d('cdata',rescale(M),'texture','2D');
    view(3);
    axis tight; % daspect([1 1 .4])
    colormap bone(256);
    alphamap('rampup');
    % alphamap(.06 .* alphamap);
    vol3d(h);
end
if not(isempty(str))
    title(str);
end

if nargout>0 
    if exist('h')
        h1 = h;
    else 
        h1 = [];
    end
end


if nargin==5 && c==a*b
    linkaxes(axlist, 'xy');
end