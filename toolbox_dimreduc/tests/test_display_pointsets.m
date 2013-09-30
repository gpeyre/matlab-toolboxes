% load and display all point sets

% '3d_cluster'
names = {'swissroll','square','scurve','puncted_sphere','twin_peaks','toroidal_helix','gaussian', '3d_cluster'};
options.sampling = 'rand';

a = 3;
b = 3;
n = 400;

clf;
for i=1:min(a*b,length(names))
    subplot(a,b,i);
    [X,col] = load_points_set( names{i}, n, options );
    plot_scattered(X,col);
    axis off;
end