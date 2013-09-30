function V = build_alpert_matrix(x,k)

% build_alpert_matrix - build 1D Alpert basis
%   adapted for the sampling provided.
%
%   V = build_alpert_matrix(pos,k);
%
%   'pos' is a 1D vector, the location of the points.
%   'k' is the number of vanishing moments (1=>Haar, 2=>linear basis ...).
%   'V' is an orthogonal matrix, each column is a basis vector.
%
%   Copyright (c) 2004 Gabriel Peyré

n = length(x);
x = x(:);

if n<k
    k=n;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first sort the sampling, and record the order !!
[x,I] = sort(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a regroupement
options.ptmax = k;
[part,B,G] = dichotomic_grouping(x',options);

P = length(G);          % number of groups

% we have got nbr.packets = 2^(J-1)
J = log2(P)+1;

si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of moments matrices
for i = 1:P
    % build moment matrix
    xs = x( si(i):(si(i+1)-1) );        % selected points
    
    kk = G(i);                  % nbr of point in this bin, ie. 2k in regular case
    % for bookkeeping we need 2*k polynmials
    [xJ,xI] = meshgrid(0:2*k-1,0:kk-1);
    Mi{i} = xs(xI+1).^xJ;    
    % for orthogonalization the matrix must be square
    [xJ,xI] = meshgrid(0:kk-1, 0:kk-1);
    M = xs(xI+1).^xJ;    
    [Ui{i},R] = qr(M);    Ui{i} = transpose(Ui{i});
end

% initialization of first matrix
U = zeros(n,n);
for i = 1:P
    selj = si(i):(si(i+1)-1);                   % selected column
    % to keep : upper part is of size n-P*k
    offs = si(i)-(i-1)*k;   % offset on row
    long = G(i)-k;          % length on row
    U( offs + (0:long-1), selj )       = Ui{i}((k+1):end, :);   % we keep the G(i)-k last
    % to retransform : lower part is of size P*k
    offs = n - P*k+1 + (i-1)*k;
    U( offs+(0:k-1), selj )   = Ui{i}(1:k, :);   % the first k doesn't have vanishing moments, keep them   
end

for j=2:J   % for each scale
    
    % at this scale, we have nj = P/2^(j-1) groups
    nj = P/2^(j-1) ; % n/(2^j*k);
    mj = nj*2*k;     % total length of the blocks
    
    % update each sub matrix
    for i = 1:nj
              
        M = zeros(2*k, 2*k);
        M(1:k,:)        = Ui{2*i-1}(1:k,:)  *   Mi{2*i-1};         % Ui^U is just k first row
        M((k+1):2*k,:)  = Ui{2*i}(1:k,:)    *   Mi{2*i};
        MMi{i} = M;
        [UUi{i},R] = qr(MMi{i});    UUi{i} = transpose(UUi{i});
        
    end
    Mi = MMi;
    Ui = UUi;
    
    % lower part of the multiplicative matrix
    UU = zeros(mj,mj);
    for i = 1:nj
        UU( k*(i-1)+(1:k), 2*k*(i-1)+(1:2*k) )          = mlow( Ui{i} );
        UU( mj/2+k*(i-1)+(1:k), 2*k*(i-1)+(1:2*k) )     = mup( Ui{i} );    
    end
    
    % multiplicative matrix
    Uj = eye(n);
    Uj( (end-mj+1):end, (end-mj+1):end ) = UU;
    
    % update vectors
    U = Uj*U;
end

V = U';

% reorder the index of the vectors !
II = reverse_permutation(I);
V = V(II,:);


% extract uper part of the matrix
function MU = mup(M)
MU = M(1:end/2, :);

% extract lower part of the matrix
function ML = mlow(M)
ML = M((end/2+1):end, :);