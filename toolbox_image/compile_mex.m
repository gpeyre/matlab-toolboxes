clear all;
fprintf('Compiling mex files ... ');

% compile mex files
mex mex/perform_adaptive_filtering.cpp

disp('done.');