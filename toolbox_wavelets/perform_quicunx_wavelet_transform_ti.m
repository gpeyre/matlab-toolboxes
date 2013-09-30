function M = perform_quicunx_wavelet_transform_ti(M,Jmin,options)

% perform_quicunx_wavelet_transform_ti - translation invariant quincunx wavelets
%
% Forward
%   MW = perform_quicunx_wavelet_transform_ti(M,Jmin,options);
% Backward
%   M = perform_quicunx_wavelet_transform_ti(MW,Jmin,options);
%
%   The implementation is based on lifting.
%
%   You can set the number of primal (analysis) and dual (synthesis)
%   vanishing moments
%       options.primal_vm (only 2/4 is supported).
%       options.dual_vm (only 2/4 is supported).
%
%   You can set the boundary conditions to 
%       options.bound = 'per'   (periodic)
%       options.bound = 'sym'   (symmetric)
%
%   This transform (but not TI) is described in
%       Wavelet Families of Increasing Order in Arbitrary Dimensions
%       Jelena Kovacevic, Wim Sweldens
%       IEEE Trans. Image Proc, 2000.
%       
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
n = size(M,1);
Jmax = log2(n)-1;
% number of scales
J = (Jmax-Jmin+1)*2;
bound = getoptions(options, 'bound', 'per');
primal_vm = getoptions(options, 'primal_vm', 4);
dual_vm = getoptions(options, 'dual_vm', 4);

[dX,dY,w] = get_quincunx_filter(primal_vm);
[dX1,dY1,w1] = get_quincunx_filter(dual_vm);

dir = +1;
if size(M,3)>1
    dir=-1;
end

if dir==1
    M0 = M;
else
    M0 = M(:,:,end);
end

[Y,X] = meshgrid(1:n,1:n);

jlist = 0:J-1;
if dir==-1
    jlist = jlist(end:-1:1);
end

W_ini = repmat( reshape(w,1,1,length(w)), n,n );
W1_ini = repmat( reshape(w1,1,1,length(w1)), n,n )/2;
Xn = W_ini; Yn = W_ini;
Xn1 = W_ini; Yn1 = W_ini;   

for j=jlist
    j1 = floor(j/2);
    dj = 2^j1;
    
    % rotate the filters
    [dXj,dYj] = rotate_quincunx(dX,dY, j);
    [dXj1,dYj1] = rotate_quincunx(dX1,dY1, j);
    
    % build set of indices 
    for k=1:length(dX)
        Xn(:,:,k) = X+dXj(k);
        Yn(:,:,k) = Y+dYj(k);
    end
    for k=1:length(dX1)
        Xn1(:,:,k) = X+dXj1(k);
        Yn1(:,:,k) = Y+dYj1(k);
    end

    % boundary conditions
    W = W_ini; W1 = W1_ini;
    if strcmp(bound,'sym')
        I = find( Xn>n | Xn<1 |Yn>n | Yn<1 ); W(I) = 0; 
        I = find( Xn1>n | Xn1<1 |Yn1>n | Yn1<1 ); W1(I) = 0;
    end
    Xn = mod(Xn-1,n)+1; Yn = mod(Yn-1,n)+1;
    Xn1 = mod(Xn1-1,n)+1; Yn1 = mod(Yn1-1,n)+1;
    
    % linear indexes
    In = Xn+(Yn-1)*n; In1 = Xn1+(Yn1-1)*n;

    if dir==1
        % detail coefficients
        d = sum(W,3); d(d==0) = 1;  % 0 should not happend anyway
        D = M0 - sum(M0(In).*W,3) ./ d;
        M(:,:,j+1) = D;
        % update coarse
        d = sum(W1,3); d(d==0) = 1;
        M0 = M0 + .5 * sum(D(In1).*W1,3) ./ d;
    else
        %% NB : since we are in TI mode, needs to be carefull and mix 2
        %% retrieves
        % retrieve coarse
        D = M(:,:,j+1);
        d = sum(W1,3); d(d==0) = 1;
        M0a = M0 - .5 * sum(D(In1).*W1,3) ./ d;
        % other retrieve
        d = sum(W,3); d(d==0) = 1;
        M0 = M(:,:,j+1) + sum(M0a(In).*W,3) ./ d;
        M0 = (M0+M0a)/2;
    end
    
    
end
if dir==1
    % record coarse
    M(:,:,end+1) = M0;
else
    M = M0;
end


%%%%%%%%%%%%%%%%%
function [dX,dY,w] = get_quincunx_filter(vm)

switch vm
    case 2
        dX = [0  0 1 -1];
        dY = [1 -1 0  0];
        w  = [1 1 1 1];
    case 4
        dX = [0  0 1 -1 2 2 -2 -2 1 -1 1 -1];
        dY = [1 -1 0  0 1 -1 1 -1 2  2 -2 -2];
        w  = [10 10 10 10 -1 -1 -1 -1 -1 -1 -1 -1];
    case 6
        dX = [0  0 1 -1 2 2 -2 -2 1 -1 1  -1   0 0 3 -3 ...
                3 3 -3 -3 2 2 -2 -2];
        dY = [1 -1 0  0 1 -1 1 -1 2  2 -2 -2   3 -3 0 0 ...
                2 -2 2 -2 3 -3 3 -3];
        w  = [174*ones(1,4) -27*ones(1,8) 2*ones(1,4) 3*ones(1,8)];
    otherwise
        error('Only 2/4 vanishing moments are supported.');
        
end
w = w/sum(w);



%%%%%%%%%%%%%%%%%
function [dXj,dYj] = rotate_quincunx(dX,dY, j)
A = [1 1;-1 1]; A = A^j;
dXj = A(1,1)*dX + A(1,2)*dY;
dYj = A(2,1)*dX + A(2,2)*dY;
