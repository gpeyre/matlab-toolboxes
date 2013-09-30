function M = compute_dead_leaves_image(n,sigma,options)

% compute_dead_leaves_image - compute a random image using the dead-leaves model
%
%   M = compute_dead_leaves_image(n,sigma,options);
%
%   n is the size of the image.
%   sigma>0 control the repartition of the size of the basic shape:
%       sigma --> 0  gives more uniform repartition of shape
%       sigma=3 gives a nearly scale invariant image.
%   options.nbr_iter put a bound on the number of iteration.
%   options.shape can be 'disk' or 'square'
%
% References : 
%   Dead leaves correct simulation :
%    http://www.warwick.ac.uk/statsdept/staff/WSK/dead.html
%
%   Mathematical analysis
%    The dead leaves model : general results and limits at small scales
%    Yann Gousseau, Fran¸cois Roueff, Preprint 2003
%
%   Scale invariance in the case sigma=3
%    Occlusion Models for Natural Images: A Statistical Study of a Scale-Invariant Dead Leaves Model
%     Lee, Mumford and Huang, IJCV 2003.
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;

if nargin<2
    sigma = 3;
end
if isfield(options,'rmin')
    rmin = options.rmin;
else
    rmin = 0.01;    % maximum proba for rmin, shoult be >0
end
if isfield(options,'rmax')
    rmax = options.rmax;
else
    rmax = 1;    % maximum proba for rmin, shoult be >0
end
if isfield(options,'nbr_iter')
    nbr_iter = options.nbr_iter;
else
    nbr_iter = 5000;
end
if isfield(options,'shape')
    shape = options.shape;
else
    shape = 'disk';
end

M = zeros(n)+Inf;

x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);


% compute radius distribution
k = 200;        % sampling rate of the distrib
r_list = linspace(rmin,rmax,k);
r_dist = 1./r_list.^sigma;
if sigma>0
    r_dist = r_dist - 1/rmax^sigma;
end
r_dist = rescale( cumsum(r_dist) ); % in 0-1

m = n^2;

for i=1:nbr_iter
    

	% compute scaling using inverse mapping
    r = rand(1);  
	[tmp,I] = min( abs(r-r_dist) );
    r = r_list(I);
    
    x = rand(1);    % position 
    y = rand(1);
    a = rand(1);    % albedo
    
    if strcmp(shape, 'disk')
        I = find(isinf(M) & (X-x).^2 + (Y-y).^2<r^2 );
    else
        I = find(isinf(M) & abs(X-x)<r & abs(Y-y)<r );
    end
    
    m = m - length(I);
    M(I) = a;
    
    if m==0
        % the image is covered
        break;
    end
end

% remove remaining background
I = find(isinf(M));
M(I) = 0;