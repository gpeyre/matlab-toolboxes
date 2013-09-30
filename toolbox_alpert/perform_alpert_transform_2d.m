function [w,info] = perform_alpert_transform_2d(v,pos,alpert_vm, dir, options)

% perform_alpert_transform_2d - transform a 2D signal.
%
%
%   [w,info] = perform_alpert_transform_2d(v,pos,alpert_vm, dir, options);
%
%   'v' is a 1D vector, the value of the function at each sampling location.
%   'pos' is a 2D vector, pos(:,i) is the ith point.
%   'alpert_vm' is the number of vanishing moments (1=>Haar, 2=>linear basis ...).
%       * 'alpert_vm' can be an integer, and then the algorithm will use the same
%         order for X and Y direction.
%       * 'alpert_vm' can be a couple of integer alpert_vm=[kx,ky] and 'kx' will be the
%         order on the X direction, and 'ky' the order on the Y direction.
%       * 'alpert_vm' can be a set of monomial, see below ('degree_type') for
%         further comments.
%   'dir' is 1 for fwd transform and -1 for bwd.
%
%   'options' is a structure that can contains the followind field :
%   'degree_type' : Polynomial degree. By default, the multiresolution spaces are defined
%       as piecewise polynimals P that satisfy
%           degX(P)<alpert_vm(1) and degY(P)<alpert_vm(2).
%       If you want to use spaces defined by 
%           degX(P)+degY(P) < alpert_vm(1)
%       then you should specify degree_type='sum'
%       (default is degree_type='max').
%       If you want to define your own multiresolution space, 
%       you can provide your own monomials in 'alpert_vm' and then 
%       set degree_type='user_defined'. It's a bit tricky
%       because you have to provide an even number of monomials, 
%       twice more than needed. Suppose you want to use
%       as multiresolution basis the polynomials {1,X}, then you can set 
%           alpert_vm = [[0;0],[1;0],[0;1],[1;1]];
%   'part_type': for automatic partition, the way the algorithm
%       will perform the grouping (can be either '1axis', '2axis' or 'kmeans', 
%       type 'help dichotomic_grouping' for more info).
%   'part': if you don't want to use automatic grouping, then 
%       you can provide a cell array that contains a binary grouping of the points
%       (same format as 'dichotomic_grouping' function).
%
%   'w' is the transformed data.
%   'info' is a struct containing the localisation information for each
%       basis Alpert vector.
%       'info.l' is the scale of the vector (0=coarse scale).
%       'info.n' is the space location of the vector.
%       'info.k' is the number of multiwavelet (in [1,...,alpert_vm(1)*alpert_vm(2)]).
%
%   WARNING: the function will try to use the mex-compiled function
%       perform_moment_transform.dll if possible, and then it
%       won't retrun 'info'. Otherwise, it will use the slower
%       function 'perform_moment_transform_slow' and 
%       'info' will be returned.
%
%   Copyright (c) 2004 Gabriel Peyré

v = v(:);

if nargin<2
    error('You must provide sampling location in pos.');
end
if nargin<3
    alpert_vm = 3;
end
if nargin<4
    dir=1;
end

options.null = 0;

if isfield(options, 'degree_type')
    degree_type = options.degree_type;
else
    if length(alpert_vm(:))<=2
        degree_type='max';
    else
        degree_type='user_defined';
    end
end

use_mex = getoptions(options, 'use_mex', 1);
part_type = getoptions(options, 'part_type', '2axis');
part = getoptions(options, 'part', []);


if length(alpert_vm)==1
    % use same order for X and Y
    alpert_vm = [alpert_vm alpert_vm];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of points
n = size(pos,2);
if n==0
    w = [];
    info.l = [];
    info.n = [];
    info.k = [];
    return;
end

switch lower(degree_type)
    case 'max'
        if n<alpert_vm(1)*alpert_vm(2)
            % special case, not enough data.
            % We want kx*ky<n subject to the same ratio, ie kx/ky=alpert_vm(1)/alpert_vm(2);
            kk = alpert_vm; % keep old one
            kx = floor( sqrt(n*alpert_vm(1)/alpert_vm(2)) );
            if kx==0
                kx=1;
            end
            ky = floor( n/kx );
            if ky==0
                ky=1;
            end
            alpert_vm = [kx ky]; 
        end
        % define the monomial basis
        [MX,MY] = meshgrid(0:2*alpert_vm(1)-1, 0:alpert_vm(2)-1);
        monomials = [MX(:)';MY(:)'];
    case 'sum'
        alpert_vm = alpert_vm(1);
        while n<alpert_vm*(alpert_vm+1)/2
            alpert_vm = alpert_vm-1;
        end
        % define the monomial basis
        s = 0; m = 0; i = 0;
        monomials = [];
        while s<alpert_vm*(alpert_vm+1)
            s = s+1;
            monomials = [monomials, [i;m-i]];
            i = i+1;
            if i>m
                m = m+1;
                i = 0;
            end
        end
    case 'user_defined'
        monomials = alpert_vm;
    otherwise
        error('Unknown degree type.');
end

if size(monomials,1)~=2
    monomials = monomials';
end
if size(monomials,1)~=2
    error('monomials must be of size 2x(2*k2).');
end

k2 = size(monomials,2)/2;     % equivalent to k^2 in Alpert paper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a regroupement
if isempty(part) && (~strcmp(part_type, '1axis') || exist('perform_moment_transform')==0 || use_mex==0 )
    clear options;
    options.ptmax = k2;
    options.part_type = part_type;
    [part,B,G] = dichotomic_grouping(pos,options);
end


if exist('perform_moment_transform')>0 && use_mex
    if strcmp(part_type, '1axis') && isempty(part)
        % will automatically compute the parition using X axis
        w = perform_moment_transform(v,pos, monomials, dir);
    else
        w = perform_moment_transform(v,pos, monomials, dir, part);
    end
    info = [];
else
    [w,info] = perform_moment_transform_slow(v,pos, monomials, dir, part);
end