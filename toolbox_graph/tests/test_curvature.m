% test for computation of curvature on meshes

name = 'mushroom';
name = 'fandisk';
name = 'horse';
name = 'nefertiti';
name = 'elephant-50kv';
name = 'bunny';
name = 'skull';
name = 'lion-head';
name = 'vase-lion';
name = 'david50kf';
name = 'david_head';
name = 'armadillo';
name = 'aphro';
name = 'gargoyle';
name = 'tyra';
name = 'screwdriver';
name = 'hand';
options.name = name;

path(path, '../toolbox_graph_data/off/');

rep = ['results/curvature/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

[vertex,face] = read_mesh(name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% local covariance analysis
options.covariance_smoothing = 15;
[C,U,D] = compute_mesh_local_covariance(vertex,face,vertex,options);

% options for display
tau = 1.2;
options.normal_scaling = 1.5;


options.normal = squeeze(U(:,2,:));
clf;
options.face_vertex_color = perform_saturation( -D(2,:)' - D(3,:)',tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
saveas(gcf, [rep name '-covariance.png'], 'png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% curvature

options.curvature_smoothing = 10;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);

options.normal = [];
clf;
options.face_vertex_color = perform_saturation(Cmax,tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
saveas(gcf, [rep name '-cmax.png'], 'png');
clf;
options.face_vertex_color = perform_saturation(Cmin,tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
saveas(gcf, [rep name '-cmin.png'], 'png');
clf;
options.face_vertex_color = perform_saturation(Cmean,tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
saveas(gcf, [rep name '-cmean.png'], 'png');
clf;
options.face_vertex_color = perform_saturation(Cgauss,tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
saveas(gcf, [rep name '-cgauss.png'], 'png');
clf;
options.face_vertex_color = perform_saturation(abs(Cmin)+abs(Cmax),tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
saveas(gcf, [rep name '-cabs.png'], 'png');
