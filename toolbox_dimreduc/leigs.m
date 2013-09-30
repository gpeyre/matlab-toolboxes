function [E,V] = leigs(X, NE, PARAM, TYPE)

% Laplacian Eigenmaps Algorithm
%
% please refer to University of Chicago
% Computer Science Technical Report TR-2002-01
% Mikhail Belkin, Partha Niyogi
% Laplacian Eigenmaps for Dimensionality Reduction and Data Representation
% Note that Laplacian, not normalized Laplacian is used here
% http://www.cs.uchicago.edu/research/publications/techreports/TR-2002-1
%
%
% Calculate the graph laplacian of the adjacency graph of data set X.
%
% [E,V] = leigs(X, NE, PARAM, TYPE)
%
%   X - (dimension x nbr_points) matrix.
%   NE - number of eigenvectors
%   TYPE - string 'nn' or string 'epsballs'
%   PARAM - integer if TYPE='nn', real number if TYPE='epsballs'
%
% Returns:
%   E - (NE x nbr_points) matrix with eigenfunctions,
%   V is (NE x NE) matrix with eigenvalues on the diagonal
%
% Author:
%
%   Mikhail Belkin
%   misha@math.uchicago.edu
%
% Modified by Gabriel Peyré

if nargin<3
    PARAM = 8;
end
if nargin<4
    TYPE = 'nn';
end

if size(X,1)<size(X,2)
    X = X';
end

L = laplacian(X, TYPE, PARAM);

% normalized Laplacian
% for i=1:size (L)
%   D(i,i) = sum(L(i,:));
%   if (D(i,i) ~= 0)
%      DD(i,i) = 1/sqrt(D(i,i));
%   else disp ('warning 0');
%      DD(i,i) = 0;
%   end
% end
% LL=DD*(D-L)*DD;


opts.tol = 1e-9;
opts.issym=1; opts.disp = 5;
% [E,V] = eigs(L,NE,'sm',opts);
options.tol = 1e-9; options.disp = 0; options.isreal = 1; options.issym = 1;

% [Y,eigenvals] = eigs(M,d+1,'SM',options);
[E,V] = eigs(L,NE+2,'SM',options);
% [E,V] = eigs(L,NE,'LM',options);

%A = DD*A;

E = E';
E = E(1:NE,:);



function L = laplacian(X, TYPE, PARAM)

% Calculate the graph laplacian of the adjacency graph of data set X.
%
% L = laplacian(X, TYPE, PARAM)
%
% X - NxK matrix. Data points are rows.
% TYPE - string 'nn' or string 'epsballs'
% PARAM - integer if TYPE='nn', real number if TYPE='epsballs'
%
% Returns: L, sparse symmetric NxN matrix
%
% Example:
%
% L = laplacian(X,'nn',6)
% L contains the Laplacian of the graph obtained from connecting
% each point of the data set to its 6 nearest neigbours.
%
%
% Author:
%
% Mikhail Belkin
% misha@math.uchicago.edu
%

disp('- Computing Laplacian Egenmaps Embedding...');

% calculate the adjacency matrix for X
if strcmp(TYPE,'nn')
    % use fast code
    options.exlude_self = 1;
    options.nn_nbr = PARAM;
    [A,tmp] = compute_nn_graph(X',options);
    I = find(isinf(A)); A(I) = 0;
else
    A = adjacency(X, TYPE, PARAM);
end

W = A;

% disassemble the sparse matrix
[A_i, A_j, A_v] = find(A);

for i = 1: size(A_i)
    % replece distances by 1
    % gaussain kernel can be used instead of 1:
    % W(A_i(i), A_j(i)) = exp(-A_v^2/t);
    W(A_i(i), A_j(i)) = 1;
end;

disp('- Computing Laplacian eigenvectors.');

D = sum(W(:,:),2);
L = spdiags(D,0,speye(size(W,1)))-W;


function d = L2_distance(a,b,df)
% L2_DISTANCE - computes Euclidean distance matrix
%
% E = L2_distance(A,B)
%
%    A - (DxM) matrix
%    B - (DxN) matrix
%    df = 1, force diagonals to be zero; 0 (default), do not force
%
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description :
%    This fully vectorized (VERY FAST!) m-file computes the
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example :
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Wed Oct 20 08:58:08 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3

% Copyright notice: You are free to modify, extend and distribute
%    this code granted that the author of the original code is
%    mentioned as the original author of the code.

% Fixed by JBT (3/18/00) to work for 1-dimensional vectors
% and to warn for imaginary numbers.  Also ensures that
% output is all real, and allows the option of forcing diagonals to
% be zero.

if (nargin < 2)
    error('Not enough input arguments');
end

if (nargin < 3)
    df = 0;    % by default, do not force 0 on the diagonal
end

if (size(a,1) ~= size(b,1))
    error('A and B should be of same dimensionality');
end

if ~(isreal(a)*isreal(b))
    disp('[Warning: running distance.m with imaginary numbers.  Results may be off.]');
end

if (size(a,1) == 1)
    a = [a; zeros(1,size(a,2))];
    b = [b; zeros(1,size(b,2))];
end

aa=sum(a.*a); bb=sum(b.*b); ab=a'*b;
d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);

% make sure result is all real
d = real(d);

% force 0 on the diagonal?
if (df==1)
    d = d.*(1-eye(size(d)));
end


function A = adjacency(X, TYPE, PARAM);

% Compute the adjacency graph of the data set X
%
% A = adjacency(X, TYPE, PARAM);
%
% X - NxK matrix. Data points are rows.
% TYPE - string 'nn' or string 'epsballs'.
% PARAM - integer if TYPE='nn', real number if TYPE='epsballs'.
%
% Returns: A, sparse symmetric NxN matrix of distances between the
% adjacent points.
%
% Example:
%
% A = adjacency(X,'nn',6)
%   A contains the adjacency matrix for the data
%   set X. For each point, the distances to 6 adjacent points are
%   stored. N
%
% Note: the adjacency relation is symmetrized, i.e. if
% point a is adjacent to point b, then point b is also considered to be
% adjacent to point a.
%
%
% Author:
%
% Mikhail Belkin
% misha@math.uchicago.edu
%

if (nargin < 3) | (strcmp(TYPE,'nn') & strcmp(TYPE,'epsballs')) | ~isreal(PARAM)

    disp(sprintf('ERROR: Too few arguments given or incorrect arguments.\n'));
    disp(sprintf('USAGE:\n A = laplacian(X, TYPE, PARAM)'));
    disp(sprintf('X - the data matrix. Data points are rows.'));
    disp(sprintf('Nearest neigbors: TYPE =''nn''    PARAM = number of nearest neigbors'));
    disp(sprintf('Epsilon balls: TYPE =''epsballs''    PARAM = redius of the ball\n'));
    return;
end

n = size(X,1);
disp (sprintf ('- X: %d points in %d dimensional space.',n,size (X,2)));

switch TYPE
    case {'nn'}
        disp(sprintf('- Creating the adjacency matrix. Nearest neighbors, N=%d.', PARAM));
    case{'eps', 'epsballs'}
        disp(sprintf('- Creating the adjacency matrix. Epsilon balls, eps=%f.', PARAM));
end;


A = sparse(n,n);
step = 100;


if (strcmp(TYPE,'nn'))
    for i1=1:step:n
        i2 = i1+step-1;
        if (i2> n)
            i2=n;
        end;
        XX= X(i1:i2,:);
        dt = L2_distance(XX',X');
        [Z,I] = sort ( dt,2);



        for i=i1:i2
            if ( mod(i, 500) ==0)
                disp(sprintf('- %d points processed.', i));
            end;
            for j=2:PARAM+1
                A(i,I(i-i1+1,j))= Z(i-i1+1,j);
                A(I(i-i1+1,j),i)= Z(i-i1+1,j);
            end;
        end


    end;

    % epsilon balls
else
    for i1=1:step:n
        i2 = i1+step-1;
        if (i2> n)
            i2=n;
        end;


        XX= X(i1:i2,:);
        dt = L2_distance(XX',X');
        [Z,I] = sort ( dt,2 );

        for i=i1:i2
            if ( mod(i, 500) ==0) disp(sprintf('- %d points processed.', i)); end;
            j=2;
            while ( (Z(i-i1+1,j) < PARAM))
                j = j+1;
                jj = I(i-i1+1,j);
                A(i,jj)= Z(i-i1+1,j);
                A(jj,i)= Z(i-i1+1,j);
            end;
        end

    end;

end;



