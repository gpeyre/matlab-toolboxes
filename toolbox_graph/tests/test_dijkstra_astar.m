% test of shortest path using Dijkstra/A*

name = 'erdos';
name = 'rand_network';
name = 'rand';
name = 'square';
name = 'sanfrancisco';
name = 'northamerica';
name = 'sanjoaquin';
name = 'california';
name = 'oldenburg';
name = 'rand_clusters_triangulation';
name = 'randn_triangulation';
name = 'nefertiti.off';
name = 'rand_triangulation';

options.proba = 0.1;
options.alpha = 0.3;
options.beta = 0.3;

n = Inf;
switch lower(name)
    case 'square'
        n = 50^2;
    case {'rand_triangulation','randn_triangulation', 'rand_clusters_triangulation'}
        n = 5000;
    case 'rand'
        n = 100;
    case 'erdos'
        n = 100;
    case 'rand_network'
        n = 200;
        
end

disp('Loading file');
[A,vertex] = load_graph(name,n,options);

n = min(n,size(A,1));
A = A(1:n,1:n);
vertex = vertex(:,1:n);
A = max(A,A');

% select two points
clf;
hold on;
gplot(A,vertex');
axis tight;
axis square;
axis off;
start_point = ginput(1)';
plot( start_point(1), start_point(2), 'ro'  );
end_point = ginput(1)';
plot( end_point(1), end_point(2), 'bx'  );
hold off;

% find nearest node
d = compute_distance_to_points(vertex,start_point);
[tmp,start_point] = min(d);
d = compute_distance_to_points(vertex,end_point);
[tmp,end_point] = min(d);

% compute heuristic
v = vertex';
p = v(end_point,:);
H = sqrt( (p(1)-v(:,1)).^2 + (p(2)-v(:,2)).^2 );

% perform front propagation without heuristic
clear options;
options.end_points = end_point;
options.H = [];
disp('Performing classical Dijkstra');
[D,S] = perform_dijkstra(A, start_point, options);

% perform front propagation with euclidean heuristic
clear options;
options.end_points = end_point;
options.H = H;
disp('Performing classical A*');
[D1,S1] = perform_dijkstra(A, start_point, options);

% extract paths
path  = perform_dijkstra_path_extraction(A,D,end_point);
path1 = perform_dijkstra_path_extraction(A,D1,end_point);

% plot the paths
options.point_size = 5;
options.graph_style = 'k';
options.far_point_style = '';
clf;
subplot(1,2,1);
plot_dijkstra(A, vertex, S, path, start_point,end_point, options );
axis square;
title('Classical Dijkstra');
subplot(1,2,2);
plot_dijkstra(A, vertex, S1, path1, start_point,end_point, options );
axis square;
title('Algorihtm A^*');