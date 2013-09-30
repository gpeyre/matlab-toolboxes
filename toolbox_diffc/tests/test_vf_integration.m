% test for vf integration (aka streamlines computations)

name = 'rot';
name = 'mixed';

% time step for integration
dt = 0.1;
m = 10;
% time for display
T_list = 0:10;

n = 80;
options.normalize = 1;
vf = load_flow(name, n, options);
vf = vf(11:end,:,:);

n = size(vf,1);
p = size(vf,2);

M = perform_vf_integration(vf, dt,T_list );

% plot some evolution
nb = 40;
[Y,X] = meshgrid( round(linspace(1,p,nb)), round(linspace(1,n,nb)) );

pos = [X(:), Y(:)];

clf;
hold on;
for i=1:size(pos,1)
    a = pos(i,1);
    b = pos(i,2);
    v = M(a,b,:);
    x = real(v(:)); x = (x-1)/n;
    y = imag(v(:)); y = (y-1)/p;
    plot(x, y, 'k');
end
plot_vf(vf);
hold off;