function [x,err_tv,err_l2, err_tgt] = perform_tv_projection(x0,tau,options)

% perform_tv_projection - perform correction of the image to that it minimizes the TV norm.
%
%   [x,err_tv,err_l2] = perform_tv_projection(x0,tau,options);
%
%   Perform the following projection on the TV ball of radius tau:
%       min_x |x-x0|_2   s.t.   |x|_TV <= tau
%
%   You can weight this projection using 
%       f = options.weight_l2
%       g = options.weight_tv
%   and the algorithm actually solves for 
%       min_x sum f(i)*|x(i)-x0(i)|^2
%   subject to
%       \sum_i g(i)*|grad_i(x)|<=1
%
%   The method use a subgradient projection, as explained in:
%       P.L. Combettes, J.C. Pesquet
%       "Image restoration subject to a total variation constraint"
%       IEEE Transactions on Image Processing 13 nç9 (2004), p. 1213-1222
%   The equations for the projection on the intersection of 2 half space is in
%       P.L. Combettes,
%       "A Block-Iterative Surrogate Constraint Splitting Method for Quadratic Signal Recovery"
%       IEEE Transactions on Signal Processing 51 nç7, (2003) p. 1771-1782
%
%   The number of iteration is options.niter.
%   You can provide an initial guss in options.x;
%
%   See also: compute_total_variation.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

if size(x0,1)==1 || size(x0,2)==1
    % 1D signal
    x0 = diff(x0); x0(end+1) = 0;
    x = perform_l1ball_projection(x0,tau);
    x = cumsum(x); x = x-mean(x)+mean(x0);
    err_tv = []; 
    err_l2 = [];
    return;
end

nbdims = 2;
if size(x0,3)>1
    nbdims = 3;
end

mask = getoptions(options, 'mask', []);
x = getoptions(options, 'x', x0);
if isempty(x)
    x = x0;
end
niter = getoptions(options, 'niter', 200);
f = getoptions(options, 'weight_l2', ones(size(x0)));
g = getoptions(options, 'weight_tv', ones(size(x0)));
tv_norm = getoptions(options,'tv_norm', 'l2');
xtgt = getoptions(options, 'xtgt', []);

% correction to account for weight
x0 = x0 .* sqrt(f);
x = x.*sqrt(f);


err_tv = [];
err_l2 = [];
err_tgt = [];
for i=1:niter
    progressbar(i,niter);
    % subgradient of the total variation
    t = grad(x./sqrt(f), options);    
    if strcmp(tv_norm, 'l2')
        t = my_perform_vf_normalization(t);
    elseif strcmp(tv_norm, 'linf')  % works only in 2D ...
        if nbdims==2
            tx = t(:,:,1); ty = t(:,:,2);
            I = find(abs(tx)>abs(ty));
            J = find(abs(tx)<=abs(ty));
            tx(I) = sign(tx(I)); ty(I) = 0;
            ty(J) = sign(ty(J)); tx(J) = 0;
            t = cat(3,tx,ty);
        else
            tx = t(:,:,:,1); ty = t(:,:,:,2); tz = t(:,:,:,3);
            I = find( abs(tx)>=abs(ty) & abs(tx)>=abs(tz)  );
            J = find( abs(ty)>abs(tx) & abs(ty)>=abs(tz) );
            K = find( abs(tz)>abs(tx) & abs(tz)>abs(ty) );
            tx(I) = sign(tx(I)); ty(I) = 0; tz(I) = 0;
            ty(J) = sign(ty(J)); tx(J) = 0; tz(J) = 0;
            tz(K) = sign(tz(K)); tx(K) = 0; ty(K) = 0;
            t = cat(4,tx,ty,tz);
        end
    elseif strcmp(tv_norm, 'l1')
        t = sign(t);
    end
    % multiply by TV weight
    if nbdims==2
        t = t.*repmat(g, [1 1 2]);
    else % 3D
        t = t.*repmat(g, [1 1 1 3]);        
    end
    t = -div( t, options );
    % correct to take into account twisted product
    t = t./sqrt(f);
    % gradient projection onto TV=tau
    tau1 = compute_total_variation(x./sqrt(f), options);
    % d = sum( t(:).^2 );
    d = dotp(t,t);
    if d>1e-9
        z = x - (tau1 - tau)*t / d;
    else
        z = x;
    end
    % correction to account for weight
    % x = x.*sqrt(f); z = z.*sqrt(f);
    % projection of x0 on the intersection
    pi = dotp(x0-x,x-z);
    mu = dotp(x0-x,x0-x);
    nu = dotp(x-z,x-z);
    rho = mu*nu-pi^2;
    xold = x;
    if rho==0 && pi>=0
        x = z;
    elseif rho>0 && pi*nu>=rho
        x = x0 + (1+pi/nu)*(z-x);
    elseif rho>0 && pi*nu<rho
        x = x + nu/rho*( pi*(x0-x)+mu*(z-x) );
    else
        error('PBM');
    end
        
    if not(isempty(mask))
        x = x.*mask;
    end
    
    % record errors
    err_tv(end+1) = tau1-tau;
    err_l2(end+1) = norm( x(:)-x0(:), 'fro' );
    if not(isempty(xtgt))
        err_tgt(end+1) = norm(x-xtgt, 'fro');
    end
    if tau1<tau
        break;
    end
    
end

% backward correction
x = x./sqrt(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vf = my_perform_vf_normalization(vf)

if size(vf,4)<=1
    d = sqrt( sum(vf.^2,3) );
else
    d = sqrt( sum(vf.^2,4) );
end
d(d<1e-9) = 1;
vf = prod_vf_sf(vf,1./d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = dotp(x,y)
d = sum(x(:).*y(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v2 = prod_vf_sf(v1,s)
if size(v1,4)<=1
    v2 = v1 .* repmat(s, [1 1 2]);
else
    v2 = v1 .* repmat(s, [1 1 1 3]);
end