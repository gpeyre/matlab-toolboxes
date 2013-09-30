function X = perform_omp(D,Y,options)

% perform_omp - perform orthogonal matching pursuit
%
%   X = perform_omp(D,Y,options);
%
%   D is the dictionary of size (n,p) of p atoms
%   Y are the m vectors to decompose of size (n,m)
%   X are the m coefficients of the decomposition of size (p,m).
%
%   Orthogonal matching pursuit is a greedy procedure to minimise
%       |x|_0   under the constraint that D*x=y.
%   (maximum sparsity).
%   You can stop the iteration after |x|_0=k  by setting options.nbr_max_atoms=k
%   You can set relative tolerance options.tol so that iterations stops when
%       |Y-D*X| < tol*|Y|
%
%	This code calls the fast mex version of Antoine Grolleau
%	when possible.
%
%   Copyright (c) 2006 Gabriel Peyre


if isfield(options, 'sparse_coder') && strcmp(options.sparse_coder, 'mp')
    X = perform_mp(D,Y,options);
    return;
end

options.null = 0;
nbr_max_atoms = getoptions(options, 'nbr_max_atoms', 50);

tol = getoptions(options, 'tol', 1e-3);
use_slow_code = getoptions(options, 'use_slow_code', 0);
verb = getoptions(options, 'verb', 1);

% P : number of signals, n dimension
[n,P]=size(Y);
% K : size of the dictionary
[n,K]=size(D);


if isfield(options, 'use_mex') && options.use_mex==1 && exist('perform_omp_mex')==3
    % use fast mex interface
    X = zeros(K,P);
    for k=1:P
        X(:,k) = perform_omp_mex(D,Y(:,k),nbr_max_atoms);
    end
    return;
end

if use_slow_code
    X = perform_omp_other(Y,D,tol,nbr_max_atoms);
    return;
end

for k=1:1:P,
    if P>20 && verb==1
        progressbar(k,P);
    end
    a=[];
    x = Y(:,k);
    r = x; % residual
    I = zeros(nbr_max_atoms,1); % index of the chosen atoms
    % norm of the original signal
    e0 = sqrt( sum( r.^2 ) );
    
    for j=1:1:nbr_max_atoms,
        proj = D'*r;
        pos = find(abs(proj)==max(abs(proj)));
        pos = pos(1);
        
        I(j) = pos;        
        a=pinv(D(:,I(1:j)))*x;
        r=x-D(:,I(1:j))*a;

        % compute error
        e = sqrt( sum( r.^2 ) );
        if e/e0 < tol
            break;
        end

    end;
    
    
    temp=zeros(K,1);
    temp(I(1:length(a)))=a;
    X(:,k)=sparse(temp);
end;
return;

