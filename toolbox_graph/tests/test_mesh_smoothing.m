% test for heat diffusion on meshes

% load mesh
name = 'mushroom';
name = 'fandisk';
name = 'lion-head';
name = 'bunny';
name = 'skull';
name = 'horse';
name = 'nefertiti';
name = 'elephant-50kv';
options.name = name;
[vertex,face] = read_off([name '.off']);

rep = 'results/mesh-smoothing/';
if not(exist(rep))
    mkdir(rep);
end

nvert = size(vertex,2);
nface = size(face,2);

laplacian_type = 'distance';
laplacian_type = 'combinatorial';
laplacian_type = 'conformal';

%% symmetric laplacian
if not(strcmp(laplacian_type,'conformal'))
    options.symmetrize = 1;
    options.normalize = 0;
    L0 = compute_mesh_laplacian(vertex,face,laplacian_type,options);
    G0 = compute_mesh_gradient(vertex,face,laplacian_type,options);
    disp(['Error (should be 0): ' num2str(norm(L0-G0'*G0, 'fro')) '.']);
    options.normalize = 1;
    L1 = compute_mesh_laplacian(vertex,face,laplacian_type,options);
    G1 = compute_mesh_gradient(vertex,face,laplacian_type,options);
    disp(['Error (should be 0): ' num2str(norm(L1-G1'*G1, 'fro')) '.']);
end

%% un-symmetric laplacian
options.symmetrize = 0;
options.normalize = 1;
disp('--> Computing laplacian');
L = compute_mesh_laplacian(vertex,face,laplacian_type,options);

%% heat diffusion flow
Tlist = [0 10 40 200];
options.dt = 0.3;
clf;
for i=1:length(Tlist)
    options.Tmax = Tlist(i);
    vertex1 = perform_mesh_heat_diffusion(vertex,face,L,options);
    % display
    subplot(1,length(Tlist),i);
    plot_mesh(vertex1,face,options);
    shading interp; camlight; axis tight;
end
saveas(gcf, [rep name '-smoothing-' laplacian_type '.png'], 'png');

%% quadratic regularization
clf;
for i=1:length(Tlist)
    At = speye(nvert) + Tlist(i)*L;
    vertex1 = (At\vertex');
    % display
    subplot(1,length(Tlist),i);
    plot_mesh(vertex1,face,options);
    shading interp; camlight; axis tight;
end
saveas(gcf, [rep name '-quadreg-' laplacian_type '.png'], 'png');

%% non linear L1 minimization
options.symmetrize = 1;
options.normalize = 0;
disp('--> Computing gradient');
G = compute_mesh_gradient(vertex,face,laplacian_type,options);
Tlist = [0 0.01 0.02 0.05 0.1 0.2];
clf;
for i=1:length(Tlist)
    t = Tlist(i);
    vertex1 = vertex;
    options.niter = 500;
    options.lambda_max = t*2;
    options.lambda_min = t;
    if t>0
    for k=1:3
        [a,tmp,E] = perform_analysis_regularization(vertex(k,:)', G, options);
        vertex1(k,:) = a(:)';
    end
    end
    % display
    subplot(2,ceil(length(Tlist)/2),i);
    plot_mesh(vertex1,face,options);
    shading interp; camlight; axis tight;
end
saveas(gcf, [rep name '-l1reg-' laplacian_type '.png'], 'png');