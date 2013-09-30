clear all;
fprintf('Compiling mex files ... ');

% Orthogonal matching pursuit
mex mex/mat_omp.c mex/perform_omp.c mex/invert.c mex/matrix_vector.c -o perform_omp_mex

disp('done.');