function Mi = perform_quincunx_interpolation(M,vm)


% perform_quincunx_interpolation - interpolate data on a quincunx grid
%
%   Mi = perform_quincunx_interpolation(M,vm);
%
%    M is a 2D image where only pixels M(i+j,i-j) are available
%   (checkboard subsampling)
%
%   vm is 2, 4 or 6
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<2
    vm = 2;
end

n = size(M,1);

%% compute mask
qsub = 2;
x = floor(0:1/(qsub-1):n); x= x(1:n);
[Y,X] = meshgrid(x,x);
mask = mod( X+Y,2 )==0;

%% compute filters
[dx,dy,w] = get_quincunx_filter(vm);


%% do averaging
[Y,X] = meshgrid(1:n,1:n);
Mi = zeros(n);
[Y,X] = meshgrid(1:n,1:n);
for i=1:length(dx)
    Xi = X+dx(i); Xi(Xi<1) = 2 - Xi(Xi<1); Xi(Xi>n) = 2*n - Xi(Xi>n);
    Yi = Y+dy(i); Yi(Yi<1) = 2 - Yi(Yi<1); Yi(Yi>n) = 2*n - Yi(Yi>n);
    Mi = Mi + w(i) * M(Xi+(Yi-1)*n);
end
Mi(mask==0) = M(mask==0);


%%
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
        error('Only 2/4/6 vanishing moments are supported.');
        
end
w = w/sum(w);