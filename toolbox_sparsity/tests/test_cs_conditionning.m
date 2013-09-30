% test for CS minoration

path(path, 'toolbox/');

% number of measurements
n = 280;
% number of atoms
m = 300;
% sparsity
s = 5;

% random measurements matrix
D = compute_compressed_sensing_matrix(n,m);

Delta = compute_cs_bounds(D,s)