% test for dictionary learning with missing data
% 
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');
name = 'barb';

% size of the image to inpaint
n0 = 200;

redun = 2;        % redundandy of the dictionary
overtraining = 30;   % over training factor
        
w = 8; n = w^2;
m = round( n*redun );           % number of atoms
p = round( overtraining*m );    % number of samples
% load random patches
M0 = load_image(name);
M0 = rescale( crop(M0,n0) );
% destroy pixels
Mmask = rand(size(M0))>.3;
M = M0.*Mmask;

%% input exemplar for learning
H = compute_all_patch(M,w);
Hmask = compute_all_patch(Mmask,w);
H = reshape(H, n,size(H,4));
Hmask = reshape(Hmask, n,size(Hmask,4));
s = std(H); [tmp,sel] = sort(s); sel = sel(end:-1:1);
Y = H(:,sel(1:4:end)); Y = Y(:,1:p);
Ymask = Hmask(:,sel(1:4:end)); Ymask = Ymask(:,1:p);
Y = Y ./ repmat( sqrt(sum(Y.^2,1)), n,1 );

%% dictionary learning
options.K = m; % number of atoms
options.sparse_coder = 'omp';
options.sparse_coder = 'mp';
options.nbr_max_atoms = 4;
options.niter_grad = 60;
options.lambda_grad = .01; % gradient descent step
options.mask = Ymask;
options.niter_learning = 100;
options.use_bootstrapping = 1;
[D,X,err] = perform_dictionary_learning(Y,options);

%% display learned dictionary
nb = [10 floor(m/10)]; ndim = 2;
A = display_dictionnary(D, X, nb,ndim );
clf;
imageplot(A);

%% perform inpainting
niter_inpainting = 15;
options.sparse_coder = 'mp';
M1 = M;
for i=1:40
    M1 = perform_blurring(M1,1.2);
    M1(Mmask) = M(Mmask);
end
for i=1:niter_inpainting
    progressbar(i,niter_inpainting);
    Y1 = compute_all_patch(M1,w);
    Y1 = reshape(Y1, n,size(Y1,4));
    X1 = perform_omp(D,Y1,options);
    Y1 = reshape(D*X1,[w w 1 size(Y1,2)]);
    M1 = compute_all_patch(Y1,size(M1,1));
    % impose boundary pixels
    M1(Mmask) = M(Mmask);
end

%% display result
rep = 'results/inpainting/';
if not(exist(rep))
    mkdir(rep);
end
MM = M; MM(Mmask==0) = 0;
imageplot({M0,MM,A,M1}, {'Ground trust', 'Input','Dictionary','Inpainted'}, 2,2,1);
saveas(gcf, [rep name '-inpainted.png'], 'png');
