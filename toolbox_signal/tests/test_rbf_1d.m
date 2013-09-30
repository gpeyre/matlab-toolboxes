% test for radial basis function interpolation in 1D.

% number of points
m = 20;

n = 2024;

scaling = 2;
x = rand(m,1)*(scaling-1)+1;
x = rescale(cumsum(x));


xi = linspace(0,1,n)';


% random values
d = rescale( rand(m,1), .1,.9 );

type = {'abs', 'poly3', 'sqrt', 'gwenn'};

clf;
for i=1:length(type)

    options.rbf = type{i};
    v = perform_rbf_interpolation(xi,x,d, options);
    subplot(length(type),1,i);
    hold on;
    plot(xi, v); plot(x, d, 'r.');
    axis([0 1 0 1]);
    hold off;
    title(['RBF ' type{i}]);

end