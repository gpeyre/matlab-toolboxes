% test for the constrained delaunay triangulation

if 0
    n = 10;
    nodes = [];
    boundary = rand(2,n);
    holes = [];
else
    % Construct the L-shaped domain.
    nodes = [ 0 0; 0 2; 1 2; 1 1; 2 1; 2 0 ];
    boundary = [ 1 2; 2 3; 3 4; 4 5; 5 6; 6 1 ];
    holes = [];
    % Add random points to the domain interior.
    n = size(nodes,1);
    nbr = 500;
    for i=1:nbr
        point = 2 * rand(1,2);
        if point(1) <= 1.0 | point(2) <= 1.0
            n = n + 1;
            nodes = [ nodes; point ];
        end
    end
end


% Triangulate the domain and plot the result.
tic;
[ nodes, triangles ] = triangulate(nodes, boundary, holes);
toc;
plot_mesh(nodes, triangles);