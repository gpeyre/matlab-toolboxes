% perform a propagation on a 2D grid
m = 20;
[A,vertex] = gen_square_graph( m );
n = size(A,1);

% add little jitter to impose uniqueness
A = A + A .* rand(n,n) * 0.01;

start_points = 0.5 * n;
end_points = 4.5*m;
options.end_points = end_points;

% use fast C code
options.use_mex = 1;
[D,S] = perform_dijkstra(A, start_points, options);
% use slow Matlab code
options.use_mex = 0;
[D1,S1] = perform_dijkstra(A, start_points, options);
% use old C code 
D2 = dijkstra_fast(A, 1:length(A));

% compare mex/non-mex
I = find( D~=Inf);
% should be 0
e = norme( D(I)-D1(I) ) + norme(S-S1)

% plot the path
path = perform_dijkstra_path_extraction(A,D,end_points);
plot_dijkstra(A, vertex, S, path, start_points,end_points, options );