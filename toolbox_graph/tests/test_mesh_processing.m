% test for several mesh processing

path(path, '../toolbox_graph_data/off/');

name = 'fandisk';
name = 'skull';
name = 'lion-head';
name = 'horse';
name = 'hammerheadtriang';
name = 'elephant-50kv';
name = 'david50kf';
name = 'bunny';
name = 'mushroom';
name = 'nefertiti';
clear options;
options.name = name;
[vertex,face] = read_off([name '.off']);

rep = ['results/mesh-processing/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

nvert = size(vertex,2);
nface = size(face,2);

zoomf = 1;
if strcmp(name, 'david50kf')
    zoomf = .8;
end

vertex = vertex - repmat(mean(vertex,2), [1 nvert]);
vertex = vertex ./ repmat(max(abs(vertex(:))), size(vertex));


laplacian_type = 'conformal';
laplacian_type = 'combinatorial';
laplacian_type = 'distance';

options.symmetrize = 1;
options.normalize = 0;
L0 = compute_mesh_laplacian(vertex,face,laplacian_type,options);
options.symmetrize = 0;
options.normalize = 1;
L = compute_mesh_laplacian(vertex,face,laplacian_type,options);



clf;
plot_mesh(vertex,face,options);
shading interp;
camlight; axis tight; zoom(zoomf);
saveas(gcf, [rep name '-mesh.png'], 'png');

zoomlist = linspace(1,4,6);
for i=1:length(zoomlist)
    clf;
    plot_mesh(vertex,face,options);
    shading faceted;
    camlight; axis tight; zoom(zoomlist(i));
    saveas(gcf, [rep name '-mesh-faceted-' num2str(i) '.png'], 'png');
end

% normal displacements
options.face_vertex_color =  [];
ndepl = 6; delta = .08/ndepl;
vertex1 = vertex;
for i=1:ndepl
    progressbar(i,ndepl);
    vertex1 = perform_normal_displacement(vertex1,face,delta);
    clf;
    plot_mesh(vertex1,face,options);
    shading interp;
    camlight; axis tight; zoom(zoomf);
    saveas(gcf, [rep name '-move-pos-' num2str(i) '.png'], 'png');
end
vertex1 = vertex;
for i=1:ndepl
    progressbar(i,ndepl);
    vertex1 = perform_normal_displacement(vertex1,face,-delta);
    clf;
    plot_mesh(vertex1,face,options);
    shading interp;
    camlight; axis tight; zoom(zoomf);
    saveas(gcf, [rep name '-move-neg-' num2str(i) '.png'], 'png');
end


%% display some eigenvectors
opts.disp = 0; 
fprintf('--> computing eigendecomposition ... ');
if nvert<3000
    [V,D] = eig(full(L0));
    nb = 100;
    if nvert<400
        nb = 30;
    end
else
    nb = 60;
    [V,D] = eigs(L0,nb,'SM',opts);
    V = real(V(:,end:-1:1));
end

fprintf('done.\n');
ilist = round(linspace(3,nb, 6));
tau=2.2; % saturation for display
for i=1:length(ilist)
    % subplot(1,length(ilist),i);
    v = real(V(:,ilist(i)));
    v = clamp( v/std(v),-tau,tau );
    options.face_vertex_color = v;
    clf;
    plot_mesh(vertex,face,options);
    shading interp; camlight; axis tight; zoom(zoomf);
    colormap jet(256);
    saveas(gcf, [rep name '-eigenvectors-' num2str(i) '.png'], 'png');
end

if size(V,2)==nvert
    % display spectrum
    f = vertex(1:3,:); f = rescale(f');
    pf = V'*f;
    clf;
    eta = 4;
    if nvert<400
        eta = 1.5;
    end
    plot(pf(1:round(nvert/eta),:));
    axis([1,nvert/eta,-1,2]);
    legend('X', 'Y', 'Z');
    saveas(gcf, [rep name '-spectrum-' num2str(i) '.png'], 'png');
    % test for mesh compression
    if strcmp(name, 'mushroom')
        num = [100 60 40 20 10 6 ]; % percents
    else
        num = [100 40 20 10 5 1 ]; % percents
    end
    for i=1:length(num)
        pv = V'*vertex';
        q = round(num(i)/100*nvert); % number of zeros
        pv(q+1:end,:) = 0;
        vertex1 = (V*pv)';
        options.face_vertex_color = [];
        clf;
        plot_mesh(vertex1,face,options);
        shading interp; camlight; axis tight; zoom(zoomf);
        saveas(gcf, [rep name '-compress-' num2string_fixeddigit(num(i),3) '.png'], 'png');
    end
end
    

funcnames = {'x', 'circles', 'circlesbw'};
for i=1:length(funcnames)
    func = funcnames{i};
    switch func
        case 'x'
            f = vertex(1,:);
        case 'quad'
            f = sum(vertex(1:2,:).^2,1);
        case 'circles'
            f = cos(2*pi*5*sqrt(sum(vertex(1:2,:).^2,1)));
            f = cos(2*pi*5*vertex(1,:));
        case 'circlesbw'
            f = double( cos(2*pi*5*sqrt(sum(vertex(1:2,:).^2,1)))>0 );  
            f = double( cos(2*pi*5*vertex(1,:))>0 );           
    end
    f = rescale( f(:) );
    clf;
    options.face_vertex_color = f;
    plot_mesh(vertex,face, options);
    shading interp; % lighting none;
    camlight; axis tight; zoom(zoomf);
    colormap jet(256);
    saveas(gcf, [rep name '-func-' func '.png'], 'png');
    if size(V,2)==nvert
        % display spectrum
        pf = V'*f;
        clf;
        plot(pf(1:end/4))
        axis([1,nvert/4,-5,5]);
        saveas(gcf, [rep name '-func-' func '-spectrum.png'], 'png');
    end
    %% display laplacian
    clf;
    v = L*f; v = clamp( v/std(v),-tau,tau );
    options.face_vertex_color = v;
    plot_mesh(vertex,face, options);
    shading interp; % lighting none; 
    camlight; axis tight; zoom(zoomf);
    colormap jet(256);
    saveas(gcf, [rep name '-func-' func '-laplacian.png'], 'png');
    %% heat diffusion flow
    Tlist = [0 10 40 200]/4;
    options.dt = 0.3;
    for k=1:length(Tlist)
        options.Tmax = Tlist(k);
        f1 = perform_mesh_heat_diffusion(f(:)',face,L,options);
        % display
        options.face_vertex_color = f1(:);
        clf;
        plot_mesh(vertex,face,options);
        shading interp; camlight; axis tight; zoom(zoomf);
        colormap jet(256);
        saveas(gcf, [rep name '-func-' func '-smoothing-' num2str(k) '.png'], 'png');
    end
end

% mesh denoising
normals = compute_normal(vertex,face)';
rho = randn(1,nvert)*.01;
vertex1 = vertex + repmat(rho,3,1).*normals;
%% heat diffusion flow
Tlist = [0 2 5 10 20 40]/4;
options.dt = 0.3;
options.face_vertex_color = [];
for k=1:length(Tlist)
    options.Tmax = Tlist(k);
    vertex2 = perform_mesh_heat_diffusion(vertex1,face,L,options);
    % display
    clf;
    plot_mesh(vertex2,face,options);
    shading interp; camlight; axis tight; zoom(zoomf);
    saveas(gcf, [rep name '-denoising-' num2str(k) '.png'], 'png');
end