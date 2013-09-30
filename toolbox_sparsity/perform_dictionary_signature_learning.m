function [Dimg,D,X,err,Z] = perform_dictionary_signature_learning(Y, k, options)

% perform_dictionary_signature_learning - compute an image signature
%
%   [Dimg,D,X,err,options.Z] = perform_dictionary_signature_learning(Y, k, w, options);
%
%   The algorithm is described in
%       M. Aharon and M. Elad, 
%       " Sparse and Redundant Modeling of Image Content Using an Image-Signature-Dictionary", 
%       Submitted.
%
%   Y is a set of exemplar, each Y(:,i)=p(:) is a w*w vector where p should 
%       be some 2D patch of size (w,w) extracted from some image(s).
%   k is the width of the dictionary signature.
%
%   The number of iterations is options.niter.
%   The sparse coder used is set in options.sparse_coder to either 'omp' or
%       'mp'.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

% width of the patches
ww = size(Y,1);
w = sqrt(ww);
% number of exemplar
m = size(Y,2);
d = k;
dd = d^2; % number of atoms in the dictionary
kk = k^2;

if m<2*dd
    warning('You should increase the number of exemplars.');
end

% option for OMP
if not(isfield(options, 'nbr_max_atoms'))
    options.nbr_max_atoms = 4;
end
if not(isfield(options, 'use_mex'))
    options.use_mex = 1;
end

Dimg = getoptions(options, 'Dimg', randn(k));
niter = getoptions(options, 'niter', 40);
use_cg = getoptions(options, 'use_cg', 1);
options.niter_max = getoptions(options, 'niter_cg', 10);
sparse_coder = getoptions(options, 'sparse_coder', 'omp');

Dimg = Dimg/sqrt( sum(Dimg(:).^2) );

if isfield(options, 'strict_sparsity')
    options.nbr_atoms_max = options.strict_sparsity;
end

% remove mean
mu = sum(Y);
Y = Y - repmat( mu/ww, [ww 1] );

% compute the dictionary extraction sparse matrix
if isfield(options, 'Z')
    Z = options.Z;
else
    Z = sparse(ww*dd,kk);
    num = 0;
    for x=1:d
        for y=1:d
            num = num+1;
            selx = x:x+w-1;
            sely = y:y+w-1;
            selx = mod(selx-1,k)+1;
            sely = mod(sely-1,k)+1;
            [Ya,Xa] = meshgrid(sely,selx);
            sel = Xa(:) + k*(Ya(:)-1);
            Z( (1:ww)' + (num-1)*ww + (sel-1)*dd*ww ) = 1;
        end
    end
end
err = []; Dimg = Dimg(:);
for i=1:niter
    disp(['--> Iteration ' num2str(i) '/' num2str(niter) '.']);  
    % extract dictionary and precompute the matrices
    D = reshape( Z*Dimg, ww, dd);
    % remove zero component and normalize
    D = D - repmat( mean(D), [ww 1] );
    D = D ./ repmat( sqrt( sum(D.^2) ), [ww 1] );
    % sparse coding
    disp('-> Sparse coding.');
    if strcmp(sparse_coder, 'omp')
        X = perform_omp(D,Y,options);
    else
        X = perform_mp(D,Y,options);
    end
    X = full(X);
    err(end+1) = norm(Y-D*X, 'fro');
    % update dictionary
    if not(use_cg)
        % pseudo inverion
        D = Y * pinv(X);
    else
        options.x = D';
        [D,err1] = perform_conjugate_gradient(X*X',X*Y',options);
        D = D';
    end
    % reconstruct
    Dimg = 1/ww * Z'*D(:);
end
Dimg = reshape(Dimg,k,k);
D = reshape( Z*Dimg(:), ww, dd);
D = D - repmat( mean(D), [ww 1] );
D = D ./ repmat( sqrt( sum(D.^2) ), [ww 1] );

% add the low pass DC
D = [ones(ww,1)/sqrt(ww) D];
X = [mu/sqrt(ww); X];