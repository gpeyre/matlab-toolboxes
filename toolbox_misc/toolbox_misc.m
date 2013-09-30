% toolbox_misc - a bunch of usefull functions.
%
%   List of functions:
%       - keep_above - keep only the coefficients above threshold T, set the rest to zero.
%       - keep_biggest - keep only the n biggest coef, set the rest to zero.
%       - l2error - compute the decreasing of the L^2 error from orthogonal coefficients.
%       - norme - just the L^2 norm of a vector/matrix
%       - plot_curve - plot a 2d curve.
%       - plot_scattered - plot 2D scattered data using triangulation.
%       - rescale - rescale data in [a,b]
%       - rev_sort_abs - sort by decreasing order of absolute value
%       - reverse - flip a vector.
%       - reverse_permutation - compute the inverse of a permuation
%       - save_image - save current graphic with given base name
%       - verbose - display a string.
%       - build_vandermonde_matrix - build the Vandermonde matrix associdated to the sampling x.
%       - clamp - clamp a value 
%       - mmax - maximum entry from a matrix
%       - poly_derivate - compute the derivative of the polynomial
%       - poly_root - find the roots of a polynomial 
%       - poly_val - evaluate a polynomial 
%       - remove_doublon - remove doublon from a vector. 
%
%   Copyright (c) 2004 Gabriel Peyré

disp('type ''help toolbox_misc''.');