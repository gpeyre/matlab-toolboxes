function [u,v,px,py] = perform_tv_hilbert_projection(f,kernel,lam,options)

% Aujol & Chambolle projection for solving TV-K regularization
%
%   [u,v,px,py] = perform_tv_hilbert_projection(f,kernel,lam,options);
%
%   Solve the image separation problem :
%       u = argmin_u |f-u|_K^2 + lam*|u|_TV
%
%   f = u+v
%   u is the cartoon part.
%   v is the texture part.
%
%   where |.|_K is an hilber space norm defined by some operator K
%       <f,g>_K = <f,K*g>
%       |f|_K^2 = <f,f>_K
%
%   f is the input image.
%   kernel is a callback function of the kernel operator, should be like
%       f = kernel(f, options);
%
%	Set options.use_gabor==1 to setup the Gabor Kernel automatically
%	and options.use_gabor==-1 to set up the L2 identity kernel.
%
% Original code by Guy Gilboa, modified by Gabriel Peyre
%
%   Copyright (c) 2006 Gabriel Peyre

[ny,nx]=size(f);

options.null = 0;
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 30;
end
if isfield(options, 'dt')
    dt = options.dt;
else
    dt = 0.01;
end
if isfield(options, 'p0x') && ~isempty(options.p0x)
    p0x = options.p0x;
else
    p0x=zeros(ny,nx); 
end
if isfield(options, 'p0y') && ~isempty(options.p0y)
    p0y = options.p0y;
else
    p0y=p0x;
end

use_gabor = 1;
if isfield(options, 'use_gabor')
    use_gabor = options.use_gabor;
end

px=p0x; py=p0y;
lami=1/lam;  % here the algorithm works with inverse of lambda

if use_gabor==1
    sig_k=1; freq=0.25; bias=0;
    n=11; %kernel size
    x=ones(n,1)*(-(n-1)/2:(n-1)/2); y=x';
    gs=exp(-((x/sig_k).^2+(y/sig_k).^2)/2);  % 2D gaussian
    gs=gs/sum(sum(gs)); % normalize
    r=sqrt(x.^2+y.^2); cs=cos(2*pi*r*freq); %radially symmetric
    iKx = cs.*gs; %% "Gabor filter" for inverse K
    iKx=iKx-mean(mean(iKx)); % zero mean filter
    cp=(n+1)/2; % center of kernel
    iKx(cp,cp) = iKx(cp,cp)-bias;
elseif use_gabor==-1
    iKx = 1;
end


%% perform iterations
if use_gabor~=0
    for i=1:niter,  %% do niterations
        progressbar(i,niter);
        %% compute projection
        km1div = filter2(iKx,div(px,py))-f*lam;
        %[Gx,Gy] = grad(km1div);
        Gx=gradx(km1div); Gy=grady(km1div);
        aG = sqrt(Gx.^2+Gy.^2); %abs(G)
        px = (px + dt*Gx)./(1+dt*aG);
        py = (py + dt*Gy)./(1+dt*aG);
    end % for i
    v = filter2(iKx,div(px,py))/lam;
else
    for i=1:niter,  %% do niterations
        progressbar(i,niter);
        % compute projection
        km1div = feval( kernel, div(px,py), options ) - f*lam;
        % compute gradient       
        Gx=gradx(km1div); Gy=grady(km1div);
%        [Gx,Gy] = grad(km1div);
        aG = sqrt(Gx.^2+Gy.^2); %abs(G)
        px = (px + dt*Gx)./(1+dt*aG);
        py = (py + dt*Gy)./(1+dt*aG);
    end % for i
    v = lami*feval( kernel, div(px,py), options );
end

u=f-v;

%% additional functions
%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient (forward difference)
function [fx,fy] = grad(P)
error('Should not be used.');
fx = P(:,[2:end end])-P;
fy = P([2:end end],:)-P;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Divergence (backward difference)
function M=div(px,py)
[m,n]=size(px);
M=zeros(m,n);
Mx=M; My=M;
Mx(2:m-1,1:n)=px(2:m-1,1:n)-px(1:m-2,1:n);
Mx(1,:)=px(1,:);
Mx(m,:)=-px(m-1,:);
My(1:m,2:n-1)=py(1:m,2:n-1)-py(1:m,1:n-2);
My(:,1)=py(:,1);
My(:,n)=-py(:,n-1);
M=Mx+My;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Laplacian (central difference)
function fl = lap(P)
[gx,gy]=grad(P);
fl=div(gx,gy);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient on X (forward difference)
function M=gradx(I)
[m,n]=size(I);
M=zeros(m,n);
M(1:m-1,1:n)=-I(1:m-1,:)+I(2:m,:);
M(m,1:n)=zeros(1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient on Y (forward difference)
function M=grady(I)
[m,n]=size(I);
M=zeros(m,n);
M(1:m,1:n-1)=-I(:,1:n-1)+I(:,2:n);
M(1:m,n)=zeros(m,1);



