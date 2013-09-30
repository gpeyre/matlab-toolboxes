% toolbox_diffc - a toolbox to perform differential calculus
%   on matrix.
%
% Here are the convention :
%   - an image is a 2D array, ie a 2D function.
%   - a vector fiels M is a 3D array, M(:,:,1) is the
%       x component of the vector field.
%   - a tensorial field T is 4D array, T(:,:,i,j)
%       with 1<=i,j<=2, is the component (i,j) of the tensorial field.
%
%   List of the functions:
%       -compute_diff - compute central derivative of a vector.
%       - compute_grad - compute the gradient of an image using central differences
%       - compute_hessian - compute the hessian tensorial field.
%       - compute_laplacian - compute the laplacian of an image.
%       - compute_operator_1 - compute a 1st order differential operator.
%       - compute_operator_2 - compute a 2nd order differential operator.
%       - compute_rigidity_tensor - compute the rigidity tensorial field
%       - plot_tf - plot a tensorial field.
%       - plot_vf - plot a vector field.
%       - perform_tensor_decomp - perform an eigendecomposition.
%   Type 'help <name of the function>' for more informations.
%
%   Copyright (c) 2004 Gabriel Peyré

disp('type ''help toolbox_diffc''.');