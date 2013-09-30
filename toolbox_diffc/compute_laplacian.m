function lap = compute_laplacian(M,options)

% compute_laplacian - compute the laplacian of an image.
%
% y = compute_laplacian(M,options);
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;

lap = compute_operator_2(M,[1,1,0],options);
