function V = build_alpert_matrix_2d(pos,k, options)

% build_alpert_matrix_2d - build 2D Alpert basis adapted to the sampling provided.
%
%   V = build_alpert_matrix_2d(pos,k);
%
%   The basis is hierarchical with respect to X direction ...
%
%   'pos' is a 2D vector, pos(:,i) is the ith point.
%   'k' is the number of vanishing moments (1=>Haar, 2=>linear basis ...).
%   'V' is an orthogonal matrix, each column is a basis vector.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    options.null = 0;
end

if isfield(options, 'part_type')
    part_type = options.part_type;
else
    part_type = '1axis';
end

if isfield(options, 'part')
    part = options.part;
else
    part = [];
end

use_scaling = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of points
n = size(pos,2);
if n==0
    V = [];
    return;
end

if n<k^2
    % special case, not enough data
    k = floor(sqrt(n));
end

k2 = k^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a regroupement
if isempty(part)
    clear options;
    options.ptmax = k2;
    options.part_type = part_type;
    [part,B,G] = dichotomic_grouping(pos(1,:),options);
end

P = length(G);          % number of groups

% we have got nbr.packets = 2^(J-1)
J = log2(P)+1;

si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of moments matrices
for i = 1:P
    
    seli = part{i};    
    % build moment matrix
    xs = pos( 1,seli );        % selected points, x component
    ys = pos( 2,seli );        % selected points, y component

    kk = length(seli);                   % nbr of point in this bin, ie. 2k^2 in regular case
    M = zeros( kk, 2*k2 );
    for ii=1:kk
        j = 0;
        % first part    : polynoms X^j1*Y^j2 for 0<=j1,j2<k
        % second part   : polynoms X^j1*Y^j2 for k<=j1<2k and 0<=j2<k
        for j1=1:2*k
            for j2=1:k
                j = j+1;
                M(ii,j) = xs(ii)^(j1-1) * ys(ii)^(j2-1);        
            end
        end
    end
        
    Mi{i} = M;    
    % orthogonalize
    [Ui{i},R] = qr(Mi{i});    Ui{i} = transpose(Ui{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of first matrix
U = zeros(n,n);
for i = 1:P
    selj = part{i};                   % selected column
    % to keep : upper part is of size n-P*k^2
    offs = si(i)-(i-1)*k2;   % offset on row
    long = G(i)-k2;          % length on row
    U( offs + (0:long-1), selj )       = Ui{i}((k2+1):end, :);   % we keep the G(i)-k^2 last
    % to retransform : lower part is of size P*k2
    offs = n - P*k2+1 + (i-1)*k2;
    if offs>0
        U( offs+(0:k2-1), selj )   = Ui{i}(1:k2, :);   % the first k^2 doesn't have vanishing moments, keep them   
    else
        warning('Pbm empty bin.');
        % bug ...
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:J   % for each scale
    
    % at this scale, we have nj = P/2^(j-1) groups
    nj = P/2^(j-1) ; % n/(2^j*k);
    mj = nj*2*k2;     % total length of the blocks
    
    % update each sub matrix
    for i = 1:nj
              
        M = zeros(2*k2, 2*k2);
        M(1:k2,:)           = Ui{2*i-1}(1:k2,:)  *   Mi{2*i-1};         % Ui^U is just k first row
        M((k2+1):2*k2,:)    = Ui{2*i}(1:k2,:)    *   Mi{2*i};
        MMi{i} = M;
        [UUi{i},R] = qr(MMi{i});    UUi{i} = transpose(UUi{i});
        
    end
    Mi = MMi;
    Ui = UUi;
    
    % lower part of the multiplicative matrix
    UU = zeros(mj,mj);
    for i = 1:nj
        UU( k2*(i-1)+(1:k2), 2*k2*(i-1)+(1:2*k2) )          = mlow( Ui{i} );
        UU( mj/2+k2*(i-1)+(1:k2), 2*k2*(i-1)+(1:2*k2) )     = mup( Ui{i} );    
    end
    
    % multiplicative matrix
    Uj = eye(n);
    Uj( (end-mj+1):end, (end-mj+1):end ) = UU;
    
    % update vectors
    U = Uj*U;
end

V = U';


% extract uper part of the matrix
function MU = mup(M)
MU = M(1:end/2, :);

% extract lower part of the matrix
function ML = mlow(M)
ML = M((end/2+1):end, :);