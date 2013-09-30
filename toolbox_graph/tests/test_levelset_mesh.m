% test for level set on meshes

path(path, '../toolbox_graph_data/off/');
if not(exist('name'))
name = 'elephant-50kv';
name = 'bunny';
name = 'david50kf';
name = 'hand';
end

options.name = name;
[vertex,face] = read_off([name '.off']);

func = 'xaxis';
func = 'eigen';

switch func
    case 'xaxis'
        F = rescale(vertex(1,:))';
        tau = linspace(0.05,.95, 10);
    case 'eigen'
        options.symmetrize = 1;
        options.normalize = 1;
        disp('--> Computing Laplacian');
        L = compute_mesh_laplacian(vertex,face,'conformal',options);
        opts.disp = 0;
        p = 150;
        disp('--> Extracting Eigenvectors');
        [U,D] = eigs(L,p, 'SM', opts);
        D = diag(abs(D)); [D,I] = sort(D); U = U(:,I);
        nb = 10; sel = round(linspace(10,p,nb));
        F = U(:, sel);
        tau = 0;
end

rep = 'results/levelsets-meshes/';
if not(exist(rep))
    mkdir(rep);
end

disp('--> Displaying');
for it=1:size(F,2)
    progressbar(it,size(F,2));
    f = F(:,it);
    options.niter_averaging = 3;
    f = perform_mesh_smoothing(face, vertex, f, options);
    % extract level set
    [v1,v2] = compute_levelset_mesh(vertex,face,f,tau,options);

    % display
    lw = 3;
    fw = perform_histogram_matching(f, linspace(0,1,length(f)));
    options.face_vertex_color = fw;
    clf;
    hold on;
    plot_mesh(vertex,face,options);
    for i=1:size(v1,2)
        h = plot3( [v1(1,i) v2(1,i)], [v1(2,i) v2(2,i)], [v1(3,i) v2(3,i)], 'k' );
        set(h, 'LineWidth', lw);
    end
    hold off;
    colormap jet(256);
    saveas(gcf, [rep name '-' func '-' num2string_fixeddigit(it,2), '.png'], 'png');
end