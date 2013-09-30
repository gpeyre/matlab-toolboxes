% test for dictionary learning
% The learning is done 
%   * for piecewise-constant 1D signals, 
%   * high contrast image patches, 
% so the learned basis should looks like a Wavelet basis.
%
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');
path(path, '../toolbox_image_data/');

test_type = 'signal';
test_type = 'image';

rep = 'results/dictionary-learning/';
if not(exist(rep))
    mkdir(rep);
end

redun = 1;        % redundandy of the dictionary
overtraining = 20;   % over training factor
        
switch test_type
    case 'signal'
        redun = 1;        % redundandy of the dictionary
        % load piecewise constant signals
        n = 80;  
        K = round( n*redun );           % number of atoms
        m = round( overtraining*K );    % number of samples

        % generate random data
        s = 3; % sparsity
        Y = zeros(n,m);
        for i=1:m
            sel = randperm(n); sel = sel(1:s);
            Y(sel,i) = randn(s,1);
        end
        Y = cumsum(Y);
    case 'image'
        w = 16; n = w^2;
        K = round( n*redun );           % number of atoms
        m = round( overtraining*K );    % number of samples
        % load random patches
        name_list = {'lena', 'boat', 'flinstones', 'goldhill', 'frog', 'mountain', 'zelda'};     % 'barb', 
        nbimg = length(name_list);
        Y = [];
        for i=1:nbimg
            M = load_image(name_list{i});
            % high pass filter to whiten a little the patches
            M = perform_blurring(M,15)-M;
            H = compute_random_patches(M,w,5*m, w);
            H = reshape(H, n,5*m);
            H = H ./ repmat( sqrt(sum(H.^2,1)), n,1 );
            s = std(H); [tmp,I] = sort(s); I = I(end:-1:1);
            Y = [Y, H(:,I(1:m))];
        end
        sel = randperm(size(Y,2)); sel = sel(1:m);
        Y = Y(:,sel);
end
Y = Y - repmat(mean(Y),[n 1]);
Y = Y ./ repmat( sqrt(sum(Y.^2,1)), n,1 );

% do the learning with several methods
options.K = K;
options.niter_learning = 60;
options.sparse_coder = 'omp';
options.sparse_coder = 'mp';
options.nbr_max_atoms = 5;
options.use_mex = 0;
options.options.init_dico = 'input';

%% run learning
disp('--> Dictionary learning.');
options.D = [];
options.X = [];
options.K = K;
options.learning_method = 'mod';
[D1,X1,E1] = perform_dictionary_learning(Y,options);

% plot some learned basis vector
if strcmp(test_type, 'image')
    nb = [10 14]; ndim = 2;
else
    nb = [4 6]; ndim = 1;
end

clf;
options.ndim = ndim;
options.normalization = 'clamp';
H = display_dictionnary(D1, X1, nb, options );

warning off;
imwrite(rescale(H), [rep 'dictionary-mod.png'], 'png');
warning on;

return;

options.learning_method = 'modortho';
options.K = n;
[D0,X0,E0] = perform_dictionary_learning(Y,options);
options.K = K;

options.learning_method = 'ksvd';
[D2,X2,E2] = perform_dictionary_learning(Y,options);
options.learning_method = 'randomized';
[D3,X3,E3] = perform_dictionary_learning(Y,options);
options.options.init_dico = 'rand';
options.learning_method = 'ksvd';
[D4,X4,E4] = perform_dictionary_learning(Y,options);

clf;
plot(1:length(E1), E1, 1:length(E1), E2, 1:length(E1), E3, 1:length(E1), E4 );
legend('MOD Ortho', 'MOD', 'KSVD', 'Mixed', 'KSVD+Input random');
axis tight;



clf;
display_dictionnary(D0, X0, nb, options );



return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the rate distortion
Dlist = {D1};

niter_coding = 30;
lambda_min = lambda/5;
lambda_max = lambda*4;
err  = zeros(niter_coding, length(Dlist)); 
bits = zeros(niter_coding, length(Dlist)); 
niter_inversion = 5;
thresh_type = 'hard';

err = zeros();
for ii=1:length(Dlist)

    D = Dlist{ii};
    [a,s,b] = svd(D);
    mu = 1/max(diag(s))^2;
    X = zeros(K,m);
    % sparse inversion
    for i=1:niter_coding
        progressbar(i,niter_coding);
        t = lambda_max - (i-1)/(niter_coding-1)*(lambda_max-lambda_min);
        for j=1:niter_inversion
            X = X + mu * D'*( Y - D*X );
            X = perform_thresholding( X, t*mu, thresh_type);
        end
        % quantize
        [Xv, Xq] = perform_quantization(X, t, 1);
        % compute error
        err(i,ii) = norm( Y-D*Xv, 'fro' ) / (n*m);
        bits(i,ii) = compute_entropy(Xq);
    end

end
plot(bits, -log(err));
xlabel('#bits'); ylabel('-log(err)'); 
legend('MOD', 'K-SVD');