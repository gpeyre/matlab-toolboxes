function plot_best_basis(x, btree, style);

% plot_best_basis - plot a best local DCT basis
%
%   plot_best_basis(x, btree);
%
%   x is the signal, btree the basis tree.
%   style is either 'frame' or 'tree'.
%
%   Copyright (c) 2006 Gabriel Peyré


str = 'r';
n = length(x);

if nargin<3
    style = 'tree';
    style = 'frame';
end
    
    
if strcmp(style, 'tree')
    
    hmin = max(x);
    hmax = max(x) + 0.3*(max(x)-min(x));

    D = log2(length(btree)+1)-1;
    stack = zeros(6,2^D+1);
    k = 1;
    stack(:,k) = [ 0 0 .5 0 .5 1. ]';
    hold on;
    plot(linspace(0,1,n),x);
    while(k > 0),
        d = stack(1,k); b = stack(2,k); x = stack(3,k);
        y = stack(4,k); w = stack(5,k); h = stack(6,k);
        % rescale h
        sc = (hmax-hmin)/D;
        y0 = (y+D)/D*(hmax-hmin)+hmin; % y0 = y;
        k=k-1;
        if(btree(node(d,b)) == 0) ,  % terminal node
            plot([(x-w/10) (x+w/10)],[y0 y0 ], str)
        else
            plot([(x-w/2) (x+w/2)],[y0 y0 ], str)
            % h = max(height(node(d,b)), maxheight/D/10 ) ;
            % h = log2(b+1);
            plot([(x-w/2) (x-w/2)],[y0 (y0-h*sc) ], str)
            k = k+1;
            stack(:,k) = [(d+1) (2*b) (x-w/2) (y-h) (w/2) h]';
            plot([(x+w/2) (x+w/2)],[y0 (y0-h*sc) ], str)
            k = k+1;
            stack(:,k) = [(d+1) (2*b+1) (x+w/2) (y-h) (w/2) h]';
        end
    end
    hold off;

elseif strcmp(style, 'frame')
    % resample if necessary to avoid too big images
    nmax = 4096;
    if n>nmax
        x = interp1(linspace(0,1,n),x,linspace(0,1,nmax));
        n = nmax;
    end
    btree = btree(:);
    hold on;
    t = linspace(0,1,n+1); t(end) = [];
    plot( t, x/max(abs(x)) );
    % plot axis subdivision
    prev = [1];
    j = 0;
    % plot left/right border
    rectangle('Position', [0 -1 1 2]);
    prec = [1];
    isalive = [1];
    while ~isempty(btree)
        cur = btree(1:2^j); btree(1:2^j) = [];
        for i=1:2^j
            if cur(i)==1 && prec(i)==1 && isalive(i)==1
                % plot a split
                p = (i-1/2)/2^j;
                plot([p p], [-1 1], str);
            end
        end
        dupl = floor(1:0.5:2^j+0.5)'; % duplicate indices
        isalive = isalive(dupl) & cur( dupl );
        prec = cur( dupl );
        j = j+1;
    end
    hold off;
    axis tight; axis off;
end