% for diffusion-wavelets orthogonal basis construction

name = 'mesh';
name = 'points';

% load data set
if strcmp(name,'points')
    
    n = 300;
    vertex = [randn(n,1) randn(n,1) randn(n,1)];
    for j=1:n
        vertex(j,:) = vertex(j,:) / norm(vertex(j,:));
    end

    Options.Delta   = 20;
    Options.NPoints = 3;
    T = MakeDiffusion(vertex, 'Gauss', Options);
    
else

    filename = 'nefertiti.off';
    [vertex,face] = read_off(filename);

    n = size(vertex,2);
    nface = size(face,2);

    disp('Computing diffusion matrix.');
    A = triangulation2adjacency(face);
    T = compute_diffusion_kernel(A);
end

% plot the eigenvalues
v = sort(abs(eigs(T, n-1)));
clf;
plot( 1:n-1,v(end:-1:1) );
axis tight;
title('Eigenvalues of the diffusion kernel');

V = compute_diffusion_geometry(T);

disp('Displaying eigenvectors.');
p = 3;
clf;
for i=1:p^2
    num = 10*(i-1)+1;
    c = V(:,num);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    subplot(p,p,i);

    if strcmp(name,'points')
        plot_scattered(vertex, V(:,num));
    else
        options.face_vertex_color = c;
        plot_mesh(vertex,face, options);
        shading interp;
    end
    str = ['Eigenvector n°', num2str(num+1)];
    title(str);
end
