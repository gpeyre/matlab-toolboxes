% test for triangulation parameterization
% (i.e. parameterization of an disk-shaped 3D model)
%
%   Copyright (c) 2005 Gabriel Peyre

path(path, '../toolbox_graph_data/off/');

name = 'mannequin';
name = 'nefertiti';
[vertex,face] = read_off([name '.off']);
A = triangulation2adjacency(face);


if strcmp(name, 'mannequin')
    vertex = -vertex;
end

rep = 'results/parameterization/';
if not(exist(rep))
    mkdir(rep);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original model
clf;
plot_mesh(vertex,face);
if strcmp(name, 'mannequin')
    view(28,10); camlight; 
end
% title('Original model');
saveas(gcf, [rep name '-mesh.png'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute 1-ring

boundary_types = {'circle','square','triangle'};
nbound = length(boundary_types);
lap_types = {'combinatorial','conformal'};
nlap = length(lap_types);
options.method = 'parameterization';

lw = 3;

kk = 0;
for l = lap_types
    k = 0;
    for b = boundary_types
        kk = kk+1;
        k = k+1;

        % parameters for the parameterization
        options.boundary = cell2mat(b);
        options.laplacian = cell2mat(l);

        k = k+1;
        str = sprintf('%s, %s', options.laplacian, options.boundary);
        disp(['Computing parameterization : ', str, '.']);
        xy = compute_parameterization(vertex,face,options);

        clf;
        h = plot_graph(A,xy,'k.-');
        set(h, 'LineWidth', lw);
        saveas(gcf, [rep name '-' options.laplacian '-' options.boundary '.png'], 'png');
    end
end
