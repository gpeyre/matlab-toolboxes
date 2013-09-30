% Test of spectral compression.
%
%   Taken from the article:
%
%   Karni Z. and Gotsman C.
%   Spectral Compression of Mesh Geometry.
%   Computer Graphics (Proceedings of SIGGRAPH), pp. 279-286, 2000. 
%   <http://www.cs.technion.ac.il/~gotsman/AmendedPubl/SpectralCompression/SpectralCompression.pdf>
%
%   Copyright (c) 2004 Gabriel Peyr?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display some example of 3D meshes
filename = 'nefertiti.off';
filename = 'mushroom.off';
[vertex,face] = read_off(filename);

nvert = size(vertex,2);
nface = size(face,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot eigenvectors
disp('Computing laplacian matrix.');
options.normalize = 0;
options.symmetrize = 1;
lap = compute_mesh_laplacian(vertex,face, 'combinatorial',options);

disp('Performing SVD.');
tic;
[U,S,V] = svd(full(lap));
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
    options.face_vertex_color = c;
    plot_mesh(vertex,face, options);
    str = ['Eigenvector nÁ', num2str(num+1)];
    title(str);
    shading interp; lighting none;
end

disp('Press any key to continue.');

pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot reconstruction
disp('Displaying reconstruction.');
% reconstruction 
p = 3;
nbr_max_keep = 8;
clf;
for i=1:p^2
    subplot(p,p,i);

    keep = 1+round(i*nbr_max_keep/p^2); % nbr de pourcent gard?
    vertex2 = (U'*vertex')';
    % set threshold
    vnorm = sum(vertex2.^2, 1);
    vnorms = sort(vnorm); vnorms = reverse(vnorms);
    thresh = vnorms( round(keep/100*nvert) );
    % remove small coefs
    mask = vnorm>thresh;
    for k=1:3
        vertex2(k,:) = vertex2(k,:).*mask;
    end
    % reconstruction
    vertex2 = (U*vertex2')';
  
    plot_mesh(vertex2,face); 
    axis tight;
    shading interp; camlight;
    str = [num2str(keep), '% of the coefficients'];
    title(str);
end