function H = compute_structure_tensor(M,sigma1,sigma2, options)

% compute_structure_tensor - compute the structure tensor
%
%   T = compute_structure_tensor(M,sigma1,sigma2);
%
%   sigma1 is pre smoothing width (in pixels).
%   sigma2 is post smoothing width (in pixels).
%
%   Follows the ideas of
%   U. Köthe, "Edge and Junction Detection with an Improved Structure Tensor", 
%   Proc. of 25th DAGM Symposium, Magdeburg 2003, 
%   Lecture Notes in Computer Science 2781, pp. 25-32, Berlin: Springer, 2003
%
%   You can set:
%       - options.sub=2 [default=2] to perform x2 interpolation of the
%           gradient to avoid aliasing.
%       - options.use_renormalization=1 [detault=0] to compute the tensor
%           using unit norm gradient vectors.
%       - options.use_anisotropic=1 [detaul=0] to perform anisotropic 
%           smoothing using a hour-glass shaped gaussian kernel.
%           This better capture edges anisotropy.
%
%   Copyright (c) 2006 Gabriel Peyré

n = size(M,1);
if nargin<2
    sigma1 = 1;
end
if nargin<3
    sigma2 = 3;
end

options.null = 0;
if isfield(options, 'sub')
    sub = options.sub;
else
    sub = 2;
end
if isfield(options, 'use_anisotropic')
    use_anisotropic = options.use_anisotropic;
else
    use_anisotropic = 0;
end
if isfield(options, 'use_renormalization')
    use_renormalization = options.use_renormalization;
else
    use_renormalization = 0;
end
if isfield(options, 'm_theta')
    m_theta = options.m_theta;
else
    m_theta = 16;
end
if isfield(options, 'sigmat')
    sigmat = options.sigmat;
else
    sigmat = 0.4; % angular smoothing
end

lambda = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre smoothing
m = min(21, round(n/2)*2+1);
h = compute_gaussian_filter([m m],sigma1/(lambda*n),[n n]);
M = perform_convolution(M,h);
% compute gradient
[gx,gy] = compute_grad(M);
% take the vector orthogonal to the gradient
tmp = gx; gx = gy; gy = -tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolates
if sub>1
    meth = 'cubic';
    [X,Y] = meshgrid(1:n,1:n);
    [XI,YI] = meshgrid(1:1/2:n+1/2,1:1/2:n+1/2); XI(:,end) = n; YI(end,:) = n;
    gx = interp2(X,Y,gx,XI,YI);
    gy = interp2(X,Y,gy,XI,YI);
end

if use_renormalization
    d = sqrt( gx.^2 + gy.^2 );
    I = find(d<eps); d(I) = 1;
    gx = gx./d; gy = gy./d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute tensors
sigma0 = 16; % initial isotropy
eta0 = 12; % intial anisotropy
T = zeros(n,n,2,2);
gx2 = gx.^2;
gy2 = gy.^2;
gxy = gx.*gy;

H = zeros(n*sub,n*sub,2,2);
H(:,:,1,1) = gx2; 
H(:,:,2,2) = gy2; 
H(:,:,1,2) = gxy; 
H(:,:,2,1) = gxy; 
[e1,e2,l1,l2] = perform_tensor_decomp(H);
l1 = l1*0 + sqrt(sigma0*eta0);
l2 = l2*0 + sqrt(sigma0/eta0);
H = perform_tensor_recomp(e1,e2,l1,l2);
gx2 = H(:,:,1,1); 
gy2 = H(:,:,2,2); 
gxy = H(:,:,1,2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth
if ~use_anisotropic
    % perform isotropic smoothing
    h = compute_gaussian_filter([m m],sub*sigma2/(lambda*n),[n n]);
    gx2 = perform_convolution(gx2,h);
    gy2 = perform_convolution(gy2,h);
    gxy = perform_convolution(gxy,h);
else
    % perform anisotropic hour-glass smoothing
    theta = linspace(0,pi,m_theta+1); theta(end) = [];
    for i=1:m_theta
        h = compute_hour_glass_filters( m*sub, sub*sigma2/(lambda*n), sigmat, n, theta(i) );
        Gx2(:,:,i) = perform_convolution(gx2,h);
        Gy2(:,:,i) = perform_convolution(gy2,h);
        Gxy(:,:,i) = perform_convolution(gxy,h);
    end
    % select correct location
    Theta = repmat( mod(atan2(gy,gx),pi), [1 1 m_theta] );
    theta = repmat( reshape(theta(:),[1 1 m_theta]), [sub*n sub*n 1] );
    [tmp,I] = min( abs(Theta-theta),[],3 );
    I = reshape( (1:(sub*n)^2)' + (I(:)-1)*(sub*n)^2, sub*n,sub*n );
    gx2 = Gx2(I);
    gy2 = Gy2(I);
    gxy = Gxy(I);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble tensor
H = zeros(n,n,2,2);
H(:,:,1,1) = gx2(1:sub:end,1:sub:end);
H(:,:,2,2) = gy2(1:sub:end,1:sub:end);
H(:,:,1,2) = gxy(1:sub:end,1:sub:end);
H(:,:,2,1) = H(:,:,1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = compute_hour_glass_filters( n, sigmar, sigmat, N, theta )

% sigmar is variance in radial direction (in pixels, so sigmar=4 is a good
%   choice)
% sigmat is variance in angular direction (in radian, so sigmar=0.4 is a good
%   choice)

if length(theta)>1
    for i=1:length(theta)
        h(:,:,i) = compute_hour_glass_filters( n, sigma, N, theta(i) );
    end
    return;
end


x = ( (0:n-1)-(n-1)/2 )/(N-1);
[Y,X] = meshgrid(x,x);

r = sqrt(X.^2 + Y.^2);
phi = atan2(Y,X);

f = exp( - r.^2 / (2*sigmar^2) );
f = f .* exp( - tan(phi-theta).^2 / (2*sigmat^2) );
f = f / sum(f(:));