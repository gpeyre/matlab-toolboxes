function Bcol = perform_segmentation_colorization(B)

% perform_segmentation_colorization - return a color image
%
%   Bcol = perform_segmentation_colorization(B);
%
%   Copyright (c) 2007 Gabriel Peyre

C = [[1; 0; 0] [0;1;0] [0;0;1] [1;1;0] [1;0;1] [0;1;1] [.3;.7;1]]';
ncol = size(C,1);
Bcol = C(cat(3,B,B+ncol,B+2*ncol));