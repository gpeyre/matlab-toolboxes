% toolbox_tensor_voting - A toolbox to perform Tensor Voting computations.
%
%   Tensor Voting is a unified framework developped by Medioni et Al.
%   to infer salient structures from images.
%   The key ingredient is the desciption of both input and output 
%   features as tensors.
%
%   The ouput is produced through a voting process that
%   will integrate the various features and infers long
%   curvilinear structures (e.g. it will 'fill gaps').
%
%   The input can be a grayscale image I (see for example 'test_ball')
%   and in this case the input field is just T(x,y)=I(x*y)*Id
%   where Id is the 2x2 identity matrix.
%   In this case the voting use only a 'ball voting field'
%   that can be computed with the function 'compute_ball_tf'.
%   The voting process can be accomplished directly on the image
%   I using the function 'perform_voting_ball'.
%
%   The input can be a complete tensor field.
%   This field can be computed from a grayscale image
%   using for example the structure tensor (see for example
%   'test_structure'). In this case the voting will use
%   a 'stick voting field'.
%   The voting process can be accomplished using 
%   the function 'perform_voting'.
%   
%   The main reference about Tensor Voting is the book 
%       A Computational Framework for Segmentation and Grouping 
%       Gérard Medioni, Mi-Suen Lee, and Chi-Keung Tang 
%       Elsevier 2000
%   But in fact I will rather recommand to read these two publications
%       W.S. Tong, C.K. Tang, and G. Medioni
%       First Order Tensor Voting, and Application to 3-D Scale Analysis
%       Proc. CVPR, pp. 175-182, 2001.
%   and
%       M.S. Lee and G. Medioni, 
%       Inferred Descriptions in Terms of Curves, Regions, and Junctions from Sparse, noisy, binary Data, 
%       Proc. International Workshop on Visual Form, pp. 350-367, 1997.
%
%   Copyright (c) 2004 Gabriel Peyré
