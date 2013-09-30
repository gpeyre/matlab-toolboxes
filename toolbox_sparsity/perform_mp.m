function X = perform_mp(D,Y,options)

% perform_mp - perform matching pursuit
%
%   X = perform_mp(D,Y,options);
%
%   D is the dictionary of size (n,p) of p atoms
%   Y are the m vectors to decompose of size (n,m)
%   X are the m coefficients of the decomposition of size (p,m).
%
%   Matching pursuit is a greedy procedure to minimise
%       |Y-D*X|  subject to |x|_0 <= s
%   where s=options.nbr_max_atoms.
%
%   You can set relative tolerance options.tol so that iterations stops when
%       |Y-D*X| < tol*|Y|
%
%   D can also be a callback function D(x,dir,options), where dir=1 to
%   compute D*x and dir=-1 to compute D'*x. In this case, the norm of the
%   columns of D should be 1/options.dico_normalization.
%
%   If you known that the coefficient to recover are real, set
%   options.isreal=1.
%
%   See also: perform_omp.
%
%   Copyright (c) 2007 Gabriel Peyre


options.null = 0;

k = getoptions(options, 'nbr_max_atoms', 10);
tol = getoptions(options, 'tol', 1e-3);
mu = getoptions(options, 'dico_normalization', 1);
isreal = getoptions(options, 'isreal', 1);

if isnumeric(D)
    m = size(D,2);
else
    xx = feval(D, Y(:,1), -1, options);
    m = size(xx,1);
end
n = size(Y,1);
p = size(Y,2);

% normalize the atoms as v/|v|^2
if isnumeric(D)
    D1 = D ./ repmat( sum(D.^2), [n 1] );
    D1t = D1';
else
    D1 = D;
end
    
e0 = norm(Y,'fro');

X = sparse(m,p);
for i=1:k
    % compute correlation
    if isnumeric(D)
        C = D1t*Y;
    else
        C = feval(D1, Y, -1, options)*mu^2;
    end
    if isreal
        C = real(C);
    end
    [v,I] = max(abs(C));
    % substract coefficients
    X1 = sparse(m,p);
    X1( I + (0:p-1)*m ) = C( I + (0:p-1)*m );   
    X = X+X1;
    if isnumeric(D)
        Y = Y - D*X1;
    else
        Y = Y - feval(D, full(X1), +1, options);
    end
    if norm(Y, 'fro')<tol*e0
        break;
    end
end