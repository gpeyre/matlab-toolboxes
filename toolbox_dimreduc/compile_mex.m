mex mex/addchv.c


return;

% compile the mex files needed for NN search
mex mex/NN/nn_prepare.cpp -Imex/ -I./
mex mex/nn_search.cpp
mex mex/range_search.cpp