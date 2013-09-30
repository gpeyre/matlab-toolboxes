% test for Lloyd algorithm

n = 150;

% initial positions
X = rand(2,n);

ndensity = 80;
x = linspace(-1,1,ndensity);
[a,b] = meshgrid(x,x);

density = 'linear';
density = 'uniform';
density = 'gaussian';

switch density
    case 'uniform'
        options.density = [];
    case 'linear'
        options.density = rescale(a,.01,1);
    case 'gaussian'
        D = exp( -(a.^2+b.^2)/.2 );
        options.density = rescale(D,.01,1);
end


rep = 'results/lloyd/';
if not(exist(rep))
    mkdir(rep);
end
repiter = ['results/lloyd/iter-' density '/'];
if not(exist(repiter))
    mkdir(repiter);
end

niter = 100;
for i=1:niter
    progressbar(i,niter);
    X = perform_lloyd_iteration(X, options);
    X1 = X; X1(:,end+1:end+4) = [1 0; 0 1; 0 0; 1 1]';
    T = delaunay(X1(1,:),X1(2,:));
    lw = 1.6;
    clf; hold on;
    h = voronoi(X1(1,:),X1(2,:));
    set(h, 'LineWidth', lw);
    options.col = 'r.-'; options.lw = lw;
    plot_graph(triangulation2adjacency(T), X1, options);
    axis([0 1 0 1]);
    hold off;
    pause(.01);
    if mod(i,2)==1 && i<50
        saveas(gcf, [repiter 'lloyd-' density '-' num2str(n) '-iter-' num2string_fixeddigit(i,2) '.png'], 'png');
    end
end

saveas(gcf, [rep 'lloyd-' density '-' num2str(n) '.png'], 'png');

