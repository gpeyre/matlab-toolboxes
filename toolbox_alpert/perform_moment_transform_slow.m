function [w,info] = perform_moment_transform_slow(v,pos, monomials, dir, part)

% perform_moment_transform_slow - perform the Alpert transform in 2D.
%
%   w = perform_moment_transform(v,pos, monomials, dir, part)
% 
% This is the heart of the algorithm, but it should never be called
% from the user, you should rather use the front-end function
% perform_alpert_transform_2d.
%
%   Copyright (c) 2004 Gabriel Peyré

% some constant and variables
k2 = size(monomials,2)/2;     % equivalent to k^2 in Alpert paper.
n = size(pos,2);

P = length(part);          % number of groups
for i=1:P
    G(i) = length(part{i});    
end
si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

% we have got nbr.packets = 2^(J-1)
J = log2(P)+1;

si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Sparse Matrix construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of moments matrices (j=1)
for i = 1:P
    seli = part{i};    
    % build moment matrix
    xs = pos( 1,seli );        % selected points, x component
    ys = pos( 2,seli );        % selected points, y component    
    kk = length(seli);                  % nbr of point in this bin, ie. 2k^2 in regular case
    
    if kk<k2 || kk>2*k2
        error('Uncorrect partitionning.');
    end
    
    M = zeros( kk, 2*k2 );
    for ii=1:kk
        M(ii,:) = xs(ii).^monomials(1,:) .* ys(ii).^monomials(2,:);
    end
%    M = M(:,1:kk);  % crop 
    Mi{i} = M;    
    % orthogonalize
    [Ui{i},R] = qr(Mi{i});    
    Ui{i} = transpose(Ui{i});
end
Uj{1} = Ui;

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
        [UUi{i},R] = qr(MMi{i});    
        UUi{i} = transpose(UUi{i});
    end
    Mi = MMi;
    Ui = UUi;    
    Uj{j} = Ui;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Sparse Matrix Multiplication
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.l = ones(n,1);     % scale 
info.n = ones(n,1);     % position on X
info.k = ones(n,1);    % multiwavelet #

if dir==1
    j_list= 1:J;
else    % transpose reverse order
    j_list= J:-1:1;
end

w = v;
for j=j_list
    
    Ui = Uj{j}; 
    if j==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % special treatment for 1st scale 
        r = w;
        for i = 1:P
            selj = part{i};                   % selected column
            % to keep : upper part is of size n-P*k^2
            offs = si(i)-(i-1)*k2;   % offset on row
            long = G(i)-k2;          % length on row
            % we keep the G(i)-k^2 last
            seli = offs + (0:long-1);
            U = Ui{i}((k2+1):end, :);
            if dir==1
                r(seli) = U * w(selj);
            else
                r(selj) = U' * w(seli);    
            end
            
            info.l(seli) = 1;
            info.n(seli) = i;
            info.k(seli) = 1:length(seli);
            
            % to retransform : lower part is of size P*k2
            offs = n - P*k2+1 + (i-1)*k2;
            if offs>0
                seli = offs+(0:k2-1);
                % the first k^2 doesn't have vanishing moments, keep them   
                U = Ui{i}(1:k2, :);
                if dir==1
                    r(seli) = U * w(selj);
                else
                    r(selj) = r(selj)  + U' * w(seli);    
                end
            else    % bug ...
                % warning('Problem empty bin.');
            end
        end
        w = r;
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % at this scale, we have nj = P/2^(j-1) groups
        nj = P/2^(j-1) ; % n/(2^j*k);
        mj = nj*2*k2;     % total length of the blocks
        
        % lower part of the multiplicative matrix
        r = w;
        offs = n-mj;
        for i = 1:nj
            selj = offs + 2*k2*(i-1)+(1:2*k2);
            seli = offs + k2*(i-1)+(1:k2);
            if dir==1
                r(seli) = mlow( Ui{i} ) * w(selj);
            else
                r(selj) = mlow( Ui{i} )' * w(seli);    
            end
            
            
            info.l(seli) = j;
            info.n(seli) = i;
            info.k(seli) = 1:k2;
            
            seli = offs + mj/2+k2*(i-1)+(1:k2);
            if dir==1
                r(seli) = mup( Ui{i} ) * w(selj);
            else
                r(selj) = r(selj) + mup( Ui{i} )' * w(seli);    
            end
            
            info.l(seli) = j;
            info.n(seli) = i;
            info.k(seli) = 1:k2;
        end
        w = r;
    end
    
end


% coarse scale
info.l(seli) = j+1;
info.l = max(info.l)-info.l;

% extract uper part of the matrix
function MU = mup(M)
MU = M(1:end/2, :);

% extract lower part of the matrix
function ML = mlow(M)
ML = M((end/2+1):end, :);