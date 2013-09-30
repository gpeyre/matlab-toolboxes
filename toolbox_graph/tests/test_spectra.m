% test for the spectra of some graph.
%
%   Copyright (c) 2004 Gabriel Peyré


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectra of a triangle
nsub = 5;
[vertex,face] = compute_base_mesh('triangle',nsub);
nface = size(face,1);
nvert = max(max(face));

% build the laplacian from adjacency matrix
disp('Computing Laplacian.');
A = triangulation2adjacency(face);
lap = compute_laplacian(A);

disp('Performing SVD.');
tic;
[U,S,V] = svd(lap);
disp( sprintf('CPU time : %.2f.', toc) );

disp('Displaying eigenvectors.');
p = 3;
clf;
for i=1:p^2
    num = 14*(i-1)+1;
    c = U(:,nvert-num);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    subplot(p,p,i);
    patch('vertices',vertex,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    str = sprintf('Eigenv. n°%d', num);
    title(str);
end

pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the 30th eigenvectors for various shapes
shape_types = {'square','square1','L1'};
nsub = [4,4,3];
num1 = 15;
num2 = 30;

nshape = length(shape_types);
clf;
for i=1:nshape
    shape_type = shape_types{i};
    
    [vertex,face] = compute_base_mesh(shape_type, nsub(i));
    nface = size(face,1);
    nvert = max(max(face));
    
    % plot surface
    subplot(4,nshape,i);
    plot_mesh(vertex,face);
    axis tight;
        title( 'Manual drawing' );
    
    disp('Computing Laplacian.');
    A = triangulation2adjacency(face);
    lap = compute_laplacian(A);
    

    disp('Performing SVD.');
    tic;
    [U,S,V] = svd(lap);
    disp( sprintf('CPU time : %.2f.', toc) );
    
    % compute laplacian embeding
    vertex1 = U(:,(nvert-2):(nvert-1));
    vertex1 = [vertex1,zeros(nvert,1)];
    
    % plot laplacian embedding
    subplot(4,nshape,i+nshape);
    plot_mesh(vertex1,face);
    axis tight;
    title( 'Laplacian-based drawing' );
        
    
    % plot eigenvector n°1
    subplot(4,nshape,i+2*nshape);
    c = U(:,nvert-num1);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    patch('vertices',vertex1,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    axis tight;
        title( sprintf('Eigenvector n°%d',num1) );


    % plot eigenvector n°2
    subplot(4,nshape,i+3*nshape);
    c = U(:,nvert-num2);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    patch('vertices',vertex1,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    axis tight;
        title( sprintf('Eigenvector n°%d',num2) );
    
end