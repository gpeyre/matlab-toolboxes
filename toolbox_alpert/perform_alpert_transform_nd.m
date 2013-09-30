function [w,info,part] = perform_alpert_transform_nd(v,pos,k,dir, options)

% perform_alpert_transform_nd - transform a nD signal using a nD Alpert basis.
%
%   [w,info] = perform_alpert_transform_nd(v,pos,k, dir, part)
%
%   'pos' is a 2D vector, pos(:,i) is the ith point.
%   'k' is the number of vanishing moments (1=>Haar, 2=>linear basis ...).
%       * 'k' can be an integer, and then the algorithm will use the same
%         order for X and Y, etc. direction
%       * 'k' can be a vector of integers k=[kx,ky,...] and 'kx' will be the
%         order on the X direction, and 'ky' the order on the Y direction, etc.
%   'dir' is 1 for fwd transform and -1 for bwd.
%
%   'w' is the transformed data.
%   'info' is a struct containing the localisation information for each
%       basis Alpert vector.
%       'info.l' is the scale of the vector (0=coarse scale).
%       'info.n' is the space location of the vector.
%       'info.k' is the number of multiwavelet (in [1,...,k(1)*k(2)]).
%       'info.s' is the number of the slice.
%
%   Copyright (c) 2004 Gabriel Peyré

v = v(:);

if nargin<2
    error('You must provide sampling location in pos.');
end
if nargin<3
    k = 3;
end
if nargin<4
    dir=1;
end
if nargin<5
    options.null = 0;
end

n = size(pos,2);    % number of points
d = size(pos,1);    % number of dimensions

if length(k)==1
    % use same order for X and Y, etc
    k = repmat(k, d, 1);
end

if n==0
    w = [];
    info.l = [];
    info.n = [];
    info.k = [];
    return;
end

if n<prod(k)
    % special case, not enough data.
    % We want prod(k)<n i.e. k^d=n
    k = floor( n^(1/d) );
    k = repmat(k, d, 1);
end
k2 = prod(k);             % equivalent to k^2 in Alpert 2D.


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a regroupement
if isempty(part)
    clear options;
    options.ptmax = k2;
    options.part_type = part_type;
    [part,B,G] = dichotomic_grouping(pos,options);
end


P = length(part);          % number of groups
for i=1:P
    G(i) = length(part{i});    
end
si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

% we have got nbr.packets = 2^(J-1)
J = log2(P)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Sparse Matrix construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of moments matrices (j=1)
for i = 1:P
    
    seli = part{i};
    if d>1
        posi = pos(:,seli);
    else
        posi = pos(seli);
    end
    kk = length(seli);  % nbr of point in this bin, ie. 2k^2 in regular case
    M = zeros( kk, 2*k2 );
    if d<5
        for ii=1:kk
            j = 0;
            %%%%%%% 1D %%%%%%
            for j1=1:2*k(1)
                if d==1
                    j = j+1;
                    M(ii,j) = posi(ii)^(j1-1);  
                else
                    %%%%%%% 2D %%%%%%
                    for j2=1:k(2)
                        if d==2
                            j = j+1;
                            M(ii,j) = posi(1,ii)^(j1-1) * posi(2,ii)^(j2-1);  
                        else
                            %%%%%%% 3D %%%%%%
                            for j3=1:k(3)
                                if d==3
                                    j = j+1;
                                    M(ii,j) = posi(1,ii)^(j1-1) * posi(2,ii)^(j2-1) * posi(3,ii)^(j3-1);
                                else
                                    %%%%%%% 4D %%%%%%
                                    for j4=1:k(4)
                                        j = j+1;
                                        M(ii,j) = posi(1,ii)^(j1-1) * posi(2,ii)^(j2-1) * posi(3,ii)^(j3-1) * posi(4,ii)^(j4-1);
                                    end                                        
                                end
                            end                            
                        end
                    end
                end
            end
        end
    else
        %%% in arbitrary dimension, it is fairly complex and slow ...
        dims = k; dims(1) = 2*dims(1);
        str = '[';
        for s=1:d
            str = [str, 'j', num2str(s), ' '] ;
        end
        str = [str, '] = ind2sub(dims, jj);'];
        % jj = (j1,j2,...,jd) with j1=0...2*k(1)-1, j2=0...k(2)-1, ..., jd=0...k(d)-1
        % M(ii,jj) = pos(1,ii)^j1 * pos(2,ii)^j2 * ... * pos(d,ii)^jd
        M = zeros( kk, 2*k2 );
        for ii=1:kk
            for jj=1:2*k2
                eval(str);
                a = 1;
                for s=1:d
                    str3 = ['a = a * posi(', num2str(s), ',ii)^(j', num2str(s), '-1);'];
                    eval(str3);
                end
                M(ii,jj) = a;
            end
        end   
    end
    
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