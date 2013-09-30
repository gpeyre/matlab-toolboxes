% test for computations of nn

name = 'swissroll';
n = 3000;
options.sampling = 'rand';
[X,col] = load_points_set( name, n, options );
nbr_nn = 50;

options.use_nntools = 0;
[D1,nn_list1] = compute_nn_distance(X,nbr_nn, options);
options.use_nntools = 1;
[D2,nn_list2] = compute_nn_distance(X,nbr_nn, options);

disp( ['Error (should be 0)=' num2str(norm(D1-D2) + norm(nn_list1-nn_list2),2) '.'] );

L = nn_list1(1,:);
col(L) = 2;

plot_scattered(X,col);