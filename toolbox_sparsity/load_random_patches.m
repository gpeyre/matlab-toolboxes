function Y = load_random_patches(Mlist,w,m,options)

% load_random_patches - extract randomly patches
%
% Y = load_random_patches(Mlist,w,m,options)
%
% Mlist is a set of images
% w is the width of the patches
% m is the number of patches
%
%   The patches are selected to have a large contrast.
%   They are a little high passed if options.high_pass=1.
%
%   Y is a (w^2,m) matrix of patches.
%   Each patch is zero mean and unit normed.
%
% Copyright (c) 2008 Gabriel Peyre

options.null = 0;
if not(iscell(Mlist))
    Mlist = {Mlist};
end
if isstr(Mlist{1})
    Mlist = load_image(Mlist);
end

high_pass = getoptions(options, 'high_pass', 1);
high_pass_width = getoptions(options, 'high_pass_width', 15);
over_selection = getoptions(options, 'over_selection', 5);

n = w^2;
nbimg = length(Mlist);
Y = [];
for i=1:nbimg
    k = floor(m/nbimg);
    if i==nbimg
        k = m-(nbimg-1)*k;
    end
    M = Mlist{i};
    % high pass filter to whiten a little the patches
    options.bound = 'per';
    M = perform_blurring(M,high_pass_width,options)-M;
    % extract patches
    H = compute_random_patches(M, w, over_selection*m, w);
    H = reshape(H, n, over_selection*m);
    H = H ./ repmat( sqrt(sum(H.^2,1)), n,1 );
    s = std(H); [tmp,I] = sort(s); I = I(end:-1:1);
    Y = [Y, H(:,I(1:m))];
end
% further shuffle the patches
sel = randperm(size(Y,2)); sel = sel(1:m);
Y = Y(:,sel);

Y = Y - repmat(mean(Y),[n 1]);
Y = Y ./ repmat( sqrt(sum(Y.^2,1)), n,1 );