function [minscale,crit,delta_list,lgd] = compute_minimum_scale(D, options)

% compute_minimum_scale - compute minimum deconvolution scale
%
%   D should be a (p,n) dictionary matrix. Each D(:,i) is an atom located
%   at some position i.
%
%   mscale(i) is a minimum scale as computed using a given criterion:
%       i=1 WERC criterion
%       i=2 ERC criterion
%       i=3 Fuchs criterion
%   one has mscale(i)>mscale(i+1) since the Fuch criterion is the finest
%   criterion (depends on sign).
%
%   options.mscale_type can be 'train' (spike train) or 'twodiracs'.
%   options.display=1 to display a graph of the criterions
%
%   Copyright (c) 2008 Gabriel Peyre


p = size(D,2); n = size(D,1);

options.null = 0;
mscale_type = getoptions(options, 'mscale_type', 'train');
delta_max = getoptions(options, 'delta_max', p/4);
disp = getoptions(options, 'display', 0);
delta_max = min(delta_max,p/2);
verb = getoptions(options, 'verb', 1);

subsampling = getoptions(options, 'subsampling', 1);

delta_list = 2:delta_max;
ntests = length(delta_list);

% record normalized gram
d = sqrt(sum(D.^2));
D1 = D ./ repmat(d, [size(D,1),1]);
G = abs(D1'*D1);


erc = [];
werc = [];
fuchs = [];
for i=1:ntests
    if verb
        progressbar(i,ntests);
    end
    delta = delta_list(i);
    % generate signal with spacing delta
    x = zeros(p,1); 
    switch mscale_type
        case 'train'
            x(1:delta:p-delta+1) = 1;
        case 'twodiracs'
            x(round(end/2)) = 1;
            x(round(end/2)+delta) = 1;
        otherwise
            error('Unknown type');
    end
    % compute criterias
    werc(i) = compute_werc_criterion(D1,x,G);
    erc(i) = 1 - compute_erc_criterion(D1,x);
    fuchs(i) = compute_fuchs_criterion(D1,x);    
end


crit = [werc; erc; fuchs];
lgd = {'werc', 'erc', 'fuchs'};
col = {'b', 'g', 'r'};
vmin = .7; vmax = 1.8;

if disp
    clf;
    hold on;
    aa = plot(delta_list*subsampling, crit); axis tight;
    set(aa, 'LineWidth', 2);
end
for i=1:3
    k = ntests;
    while crit(i,k)<1 && k>0
        k = k-1;
    end
    minscale(i) = delta_list(min(k+1,end))*subsampling;
    if disp
        aa = plot([1 1]*minscale(i), [vmin vmax], [col{i}(1) ':']);
        set(aa, 'LineWidth', 2);
    end
end
if disp
    axis([0 max(delta_list)*subsampling vmin vmax]);
    hold off; box on;
    legend(lgd);
end