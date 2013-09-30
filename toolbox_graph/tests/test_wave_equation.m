name = 'david50kf';
name = 'hand';
name = 'elephant-50kv';
name = 'bunny';
options.name = name;

path(path, '../toolbox_graph_data/off/');

rep = ['results/curvature/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

[vertex,face] = read_mesh(name);
n = size(vertex,2);

options.symmetrize = 0;
options.normalize = 1;
type = 'combinatorial';
type = 'conformal';
if not(exist('L'))
    L = compute_mesh_laplacian(vertex,face,type,options);
end

% initial conditions
sigma = .015;
if strcmp(name, 'bunny')
    sigma = .002;
end
npoints = 10; m =0;
f = zeros(n,1);
while m<npoints
    clf;
    options.face_vertex_color = f;
    plot_mesh(vertex,face,options);
    colormap jet(256);
    disp('Click on mesh and then hit enter.');
    pause; p = select3d;
    if not(isempty(p))
        d = compute_distance_to_points(vertex,p)';        
        f = f + (-1)^mod(m,2) *exp( -d/(2*sigma.^2) );
        m = m+1;
    end
end

Tmax = 80;
dt = .5;

rep = 'results/wave-equation/';
if not(exist(rep))
    mkdir(rep);
end

fprev = f;
m = 1;
for i=1:round(Tmax/dt)
    f1 = f;
    f = 2*f - fprev - dt^2 * L*f;
    fprev = f1;
    if mod(i,round(4/dt))==1
        a = f; a(a==max(a)) = max(abs(f));
        a(a==min(a)) = -max(abs(f));
        options.face_vertex_color = a;
        clf;
        plot_mesh(vertex,face,options);
        colormap jet(256);
        saveas(gcf, [rep name '-wave-eq-' num2string_fixeddigit(m,2) '.png'], 'png');
        m = m+1;
    end
end