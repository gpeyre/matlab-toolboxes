% test for the signal processing toolbox
%   
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%% 1D tests %%%%%%%%%%%%%%%%%%%%%%%
n = 128;
x = [zeros(n/2,1); ones(n/2,1)];

% different filter you can use
h = compute_gaussian_filter(13,0.1,n);
h = [0.2,0.6,0.2];
h = [0.1,0.2,0.4,0.2,0.1];
% more difficult to handle : even length works also, but it shift results by
% 1/2
h = [0.1,0.4,0.4,0.2];

y = perform_convolution(x,h);

%%%%%%%%%%%%%%%%%%%%% 2D tests %%%%%%%%%%%%%%%%%%%%%%%
n = 128;
x = -1:2/(n-1):1;
[Y,X] = meshgrid(x,x);
M = (2*(Y>=0)-1).*(2*(X>=0)-1);
M = double(M);

% odd
h = [[0 1 0];[1 2 1];[0 1 0]];
% even
h = [[0 1 1 0];[1 2 2 1];[1 2 2 1];[0 1 1 0]];
h = h/sum(sum(h));

y = perform_convolution(M,h);