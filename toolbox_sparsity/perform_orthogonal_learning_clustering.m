function [D,E] = perform_orthogonal_learning_clustering(Y, options)

% perform_orthogonal_learning_clustering - learn a set of orthogonal bases
%
%   D = perform_orthogonal_learning_clustering(Y, options);
%
%   Each Y(:,i) is a exemplar vector (e.g. a patch in an image).
%
%   D(:,:,k) is the dictionary number k.
%   E(j) is the energy of clustering at iteration j.
%
%   The total number in the dictionary is options.ndico or you can provide
%   options.num where num(i) is the initial number of the dictionary for
%   exemplar i.
%
%   You can set options.learning_mode = 'sparse' for sparse dictionary
%   learning or options.learning_mode = 'pca' for simple PCA.
%
%   options.niter is the number of iteration of the clustering algorithm.
%   options.niter_learning is the number of iteration for the learning of
%       each dictionary at each iteration of the clustering.
%
%   set options.verb=0 if you do not want output during the iterations.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;

if isfield(options, 'num')
    num = options.num;
    q = max(num);
else
    q = getoptions(options, 'ndico', 1, 1);
    % initialize the sets at random
    num = floor(linspace(1,q+1,m+1)); num(end) = [];
    num = num(randperm(m));
end

learning_mode = getoptions(options, 'learning_mode', 'sparse');
options.learning_method = 'modortho';
verb = getoptions(options, 'verb', 1);

d = size(Y,1);

niter = getoptions(options, 'niter', 20);

D = []; X = []; E = [];
for i=1:niter
    if verb
        clf;
        hist(num, 1:q); axis tight;
        title(['Number of elements in each class, iteration ' num2str(i) '/' num2str(niter)]);
        drawnow;
    end
    % first learn the dictionary
    for k=1:q
        if strcmp(learning_mode, 'sparse')
            if i>1
                options.D = D(:,:,k);
            end
            D(:,:,k) = perform_dictionary_learning(Y(:,num==k),options);
        else
            D = pca(Y,d);
        end
        % coefficients
        X(:,:,k) = D(:,:,k)' * Y;
    end
    % compute the newt clustering by minimizing L1 norm
    e = sum( abs(X), 1 );
    [e,num] = min(e,[],3);
    E(end+1) = sum(e,2);
end