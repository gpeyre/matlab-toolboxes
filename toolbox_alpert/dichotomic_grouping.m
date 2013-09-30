function [part,B,G] = dichotomic_grouping(pos,options)

% dichotomic_grouping - regroup points into 2^s (s is unknown) bins of nearly equal size. 
%
%   [part,B,G] = dichotomic_grouping_1d(pos,options);
%
%   'part' is a cell array containing the membership of each group.
%   'B' is a vector containing the membership of each point.
%   'G' is the number of points in each bin
%
%   'options' is a structure that can contains : 
%       - options.part_type: the way the dyadic split is performed, can be
%       either:
%           '1axis','2axis','3axis',etc : use the first 1,2,3,etc axis 
%              to do alternatively the split. Typically, if you have
%              3D data and you want approximatively isotropic grouping, 
%              use '3axis'.
%           'kmeans' : use a modified k-means method to do the grouping.
%       - options.ptmax: maximum number of point in each final cell.
%       - options.depthmax: depth of subdivision.
%       - options.part: an initial partition to refine.
%       - options.nb_iter : nbr of iterations for the K-means like clustering.
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<2
    options.null = 0;
end

if isfield(options, 'ptmax')
    ptmax = options.ptmax;
else
    ptmax = 0;
end

if isfield(options, 'depthmax')
    depthmax = options.depthmax;
else
    if ~isfield(options, 'ptmax')
        depthmax = 5;
    else
        depthmax = Inf;
    end
end

if isfield(options, 'nb_iter')
    nb_iter = options.nb_iter;
else
    nb_iter = 10;
end

n = size(pos, 2);
d = size(pos, 1);   % number of dimensions

if isfield(options, 'part')
    part = options.part;
    cpt = floor(log2(length(part)));
else
    part = {1:n};
    cpt = 0;
end

if isfield(options, 'part_type')
    type = options.part_type;
else
    type = '1axis';
end

if length(type)>=5 && strcmp(type(2:end), 'axis')
    % retrieve the number of axis
    type = str2num(type(1));
    type = max(min(type,d), 0);     % should be in [0,d]
else
    % otherwise use kmeans grouping
    type = 0;
end

while cpt<depthmax  % to avoid infinite loop ...
    cpt = cpt+1;
    
    part_prev = part;
    
    % perform binary split
    part1 = part;
    part = {};
    for k=1:length(part1)
        P = part1{k};
        pp = floor(length(P)/2);    % half # points
        % compute sub-partition
        if type>0
            % axis-split clustering
            s = mod(cpt-1,type)+1;
            [tmp,I] = sort( pos(s,P) );
            ppart1 = I(1:pp);
            ppart2 = I(pp+1:end);
        else
            % k-mean clustering
            [ppart1,ppart2] = dist_part( pos(:,P), nb_iter );
        end
        % assign new partition
        part{2*k-1} = P(ppart1);
        part{2*k-0} = P(ppart2);
    end
    
    for i=1:length(part)
        if length(part{i})<ptmax
            part = part_prev;
            B = zeros(n,1);
            for i=1:length(part)
                B(part{i}) = i;
                G(i) = length(part{i});  
            end
            return;
        end
    end
end

B = zeros(n,1);
for i=1:length(part)
    B(part{i}) = i;
    G(i) = length(part{i});    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,J] = dist_part(X,nb_iter)

% dist_part - partition a set of points into too components.
%
%   [I,J] = dist_part(X);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    nb_iter = 5;
end

n = size(X,2);
n1 = floor(n/2);
d = size(X,1);

% try to reduce a bit the number of iterations
nb_iter = min(nb_iter,ceil(n/2-1));

% random initialization
seed_num = floor(rand(2,1)*n)+1;
if seed_num(2)==seed_num(1)
    seed_num(2) = seed_num(1) + 1;
    if seed_num(2)>n
        seed_num(2) = 1;
    end
end
seeds = X(:,seed_num);

% compute region of influence
D = sqrt( compute_distance_to_points(X,seeds) );
[tmp,A] = sort( D(1,:)-D(2,:) );
part{1}= A(1:n1); part{2} = A(n1+1:end);

[tmp,B] = min(D);

for i=1:nb_iter
    % compute region center
    for k=1:2
        % geometric barycenter
        if ~isempty(part{k})
            seeds(:,k) = sum( X(:,part{k}), 2)/length( part{k} );
        else
            warning('Empty cluster created.');
        end
    end
    % compute region of influence
    D = sqrt( compute_distance_to_points(X,seeds) );
    [tmp,A] = sort( D(1,:)-D(2,:) );
    part{1}= A(1:n1); part{2} = A(n1+1:end);
end

I = part{1};
J = part{2};