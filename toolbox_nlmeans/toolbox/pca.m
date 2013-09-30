function [Y,X1,v,Psi] = pca(X,numvecs, options)

% pca - compute the principal component analysis.
%
%   [Y,X1,v,Psi] = pca(X,numvecs)
%
%   X is a matrix of size dim x p of data points.
%   X1 is the matrix of size numvecs x p (projection on the numvect first eigenvectors)
%	Y the matrix  of size dim x numvecs of numvecs first eigenvector of the correlation matrix X*X'
%		(this matrix is computed using the traditional flipping trick if p is large).
%	v is the vector of size numvecs of eigenvalues.
%   Psi is the mean.
%
%   Warning: the mean of X is substracted before computing the covariance
%   matrix.
%
%   You can use an iterative algorithm based 
%   on expectation maximization by setting
%       options.use_em = 1;
%   if you want a fast estimation of a few eigenvectors.
%   This algorithm use the code of Sam Roweis
%       Sam Roweis, "EM Algorithms for PCA and SPCA",
%       Neural Information Processing Systems 10 (NIPS'97) pp.626-632
%       http://www.cs.toronto.edu/~roweis/code.html
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;

% compute mean
Psi = mean(X')';

if isfield(options, 'use_em') && options.use_em==1
    if isfield(options, 'iter')
        iter = options.iter;
    else
        iter = 20;
    end
    [Y,v] = empca(X,numvecs,iter);
    X1 = Y' * X;
    return;
end


dim = size(X,1);
p = size(X,2);

% substract mean
for i = 1:p
    X(:,i) = X(:,i) - Psi;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First case: few eigenvectors asked
if numvecs<=dim && dim<p
    % do not use trick
    L = X*X'; %  should be small
    [Y,v] = eig(L);
    % sort the eigenvalues
    [Y,v] = sortem(Y,v);
    Y = Y(:,1:numvecs);
    X1 = Y' * X;
    v = diag(v);
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General case

% covariance matrix
L = X'*X;
% Eigen decomposition
[Y,v] = eig(L);
% Sort the vectors/values according to size of eigenvalue
[Y,v] = sortem(Y,v);

% Convert the eigenvectors of X'*X into eigenvectors of X*X'
Y = X*Y;

% Get the eigenvalues out of the diagonal matrix and
% normalize them so the evalues are specifically for cov(X'), not X*X'.
v = diag(v);
v = v / (p-1);

% Normalize Y to unit length, kill vectors corr. to tiny evalues
num_good = 0;
for i = 1:p
    Y(:,i) = Y(:,i)/norm(Y(:,i));
    if v(i) < 0.00001
        % Set the vector to the 0 vector; set the value to 0.
        v(i) = 0;
        Y(:,i) = zeros(size(Y,1),1);
    else
        num_good = num_good + 1;
    end;
end;
if (numvecs > num_good)
    fprintf(1,'Warning: numvecs is %d; only %d exist.\n',numvecs,num_good);
    numvecs = num_good;
end;
Y = Y(:,1:numvecs);
% perform projection
X1 = (Y')*X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vectors values] = sortem(vectors, values)

%this error message is directly from Matthew Dailey's sortem.m
if nargin ~= 2
    error('Must specify vector matrix and diag value matrix')
end;

vals = max(values); %create a row vector containing only the eigenvalues
[svals inds] = sort(vals,'descend'); %sort the row vector and get the indicies
vectors = vectors(:,inds); %sort the vectors according to the indicies from sort
values = max(values(:,inds)); %sort the eigenvalues according to the indicies from sort
values = diag(values); %place the values into a diagonal matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [evec,eval] = empca(data,k,iter,Cinit)
%[evec,eval] = empca(data,k,iter,Cinit)
%
% EMPCA
%
% finds the first k principal components of a dataset 
% and their associated eigenvales using the EM-PCA algorithm
%
% Inputs:  data is a matrix holding the input data
%               each COLUMN of data is one data vector
%               NB: mean will be subtracted and discarded
%          k    is # of principal components to find
%
% optional:
%          iters is the number of iterations of EM to run (default 20)
%          Cinit is the initial (current) guess for C (default random)
%
% Outputs:  evec holds the eigenvectors (one per column)
%           eval holds the eigenvalues
%


[d,N]  = size(data);
data = data - mean(data,2)*ones(1,N);

if(nargin<4) Cinit=[]; end
if(nargin<3) iter=20; end

[evec,eval] = empca_orth(data,empca_iter(data,Cinit,k,iter));



function [C] = empca_iter(data,Cinit,k,iter)
%[C] = empca_iter(data,Cinit,k,iter)
%
% EMPCA_ITER
%
% (re)fits the model 
%
%    data = Cx + gaussian noise 
%
% with EM using x of dimension k
%
% Inputs:  data is a matrix holding the input data
%               each COLUMN of data is one data vector
%               NB: DATA SHOULD BE ZERO MEAN!
%          k    is dimension of latent variable space 
%               (# of principal components)
%          Cinit is the initial (current) guess for C
%          iters is the number of iterations of EM to run
%
% Outputs: C is a (re)estimate of the matrix C 
%             whose columns span the principal subspace 
%

% check sizes and stuff
[p,N] = size(data);
assert(k<=p);
if(isempty(Cinit))
    C = rand(p,k);
else
    assert(k==size(Cinit,2));
    assert(p==size(Cinit,1));
    C = Cinit;
end

% business part of the code -- looks just like the math!
for i=1:iter
    % e step -- estimate unknown x by random projection
    x = inv(C'*C)*C'*data;
    % m step -- maximize likelihood wrt C given these x values
    C = data*x'*inv(x*x');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evec,eval] = empca_orth(data,C)
%[evec,eval] = empca_orth(data,Cfinal)
%
% EMPCA_ORTH
%
% Finds eigenvectors and eigenvalues given a matrix C whose columns span the
% principal subspace.
%
% Inputs:  data is a matrix holding the input data
%               each COLUMN of data is one data vector
%               NB: DATA SHOULD BE ZERO MEAN!
%          Cfinal is the final C matrix from empca.m
%
% Outputs: evec,eval are the eigenvectors and eigenvalues found
%          by projecting the data into C's column space and finding and
%          ordered orthogonal basis using a vanilla pca method
%

C = orth(C);
[xevec,eval] = truepca(C'*data);
evec = C*xevec;

function [] = assert(condition,message)

if nargin == 1,message = '';end
if isempty(message),message = 'Assert Failure.'; end
if(~condition) fprintf(1,'!!! %s !!!\n',message); end



function [evects,evals] = truepca(dataset)
% [evects,evals] = truepca(dataset)
%
% USUAL WAY TO DO PCA -- find sample covariance and diagonalize
%
% input: dataset 
% note, in dataset, each COLUMN is a datapoint
% the data mean will be subtracted and discarded
% 
% output: evects holds the eigenvectors, one per column
%         evals holds the corresponding eigenvalues
%

[d,N]  = size(dataset);

mm = mean(dataset')';
dataset = dataset - mm*ones(1,N);

cc = cov(dataset',1);
[cvv,cdd] = eig(cc);
[zz,ii] = sort(diag(cdd));
ii = flipud(ii);
evects = cvv(:,ii);
cdd = diag(cdd);
evals = cdd(ii);