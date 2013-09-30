% test for dictionary signature computation

path(path, 'toolbox/');

w = 8; ww = w^2;
k = 25; % width of the signature
overtraining = 3;
m = round( overtraining*k^2 );    % number of samples
% load random patches
name_list = {'barb'}; % {'lena', 'boat', 'flinstones', 'goldhill', 'peppers'};     % 'barb',
nbimg = length(name_list);
Y = [];
for i=1:nbimg
    M = load_image(name_list{i}); M = crop(M,256);
    % high pass filter to whiten a little the patches
    M = perform_blurring(M,15)-M;
    H = compute_random_patches(M,w,5*m, w);
    H = reshape(H, ww,5*m);
    H = H ./ repmat( sqrt(sum(H.^2,1)), ww,1 );
    s = std(H); [tmp,I] = sort(s); I = I(end:-1:1);
    Y = [Y, H(:,I(1:m))];
end
sel = randperm(size(Y,2)); sel = sel(1:m);
Y = Y(:,sel);

% options for signature learning
clear options;
options.niter = 4;
options.use_mex = 0;
options.sub = 1;
options.sparse_coder = 'mp';

% perform signature learning
err = [];
for i=1:5
    [Dimg,D,X,err1] = perform_dictionary_signature_learning(Y, k, options);
    err = [err; err1(:)];
    options.Dimg = Dimg;
end