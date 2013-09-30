% test for surface compression using fixed basis.
% This is for speeding up the spectral mesh compression 
% approach (no need to extract eigenvector). 
%   
% This idea first appears in the paper
%   Karni Z. and Gotsman C.
%   3D Mesh Compression Using Fixed Spectral Bases.
%   Proceedings of Graphics Interface, Ottawa, June 2001.
% We use here however an other approach. 
% Instead of using a spectral bassis of some kwnon regular mesh, we use 
% a known basis (cosine or wavelets) of the regular grid, 
% sample it, and re-orthogonalize it.
%
%   Copyright (c) 2005 Gabriel Peyré

% sample a basis at location xy
basis_type = 'wavelet';
basis_type = 'cosine';
basis_type = 'spectral';
basis_type = 'alpert';

rep = 'off/';
name = 'mannequin.off';
name = 'nefertiti.off';
filename = [rep name];
[vertex,face] = read_off(filename);

n = size(vertex,2);

% compute position on a square [0,1]^2
if ~strcmp(basis_type, 'spectral')
    disp('Computing 2D embedding.');
    xy = compute_parametrization(vertex,face,'combinatorial','square');
end


A = zeros(n,n);
x = xy(:,1);
y = xy(:,2);
disp('Computing transform basis.');
switch lower(basis_type)
    case 'spectral'
        % this is *not* a fixed basis
        A = triangulation2adjacency(face);
        lap = compute_laplacian(A);
        disp('Performing SVD.');
        [A,S,V] = svd(lap);
        % be sure that low frequencies comes first
        U = A(:,end:-1:1);
    case 'cosine'
        % DCT-like basis, for i,j=0...n-1
        % phi_{i,j}(x,y) = cos( pi*(i+1/2)*x ) * cos( pi*(j+1/2)*x )
        k = 0;
        for s=0:N-1 % sum of both frequencies
            if k>n
                break;
            end 
            for j=0:s
                if k>n
                    break;
                end
                k = k+1;
                i = s-j;
                A(:,k) = cos( pi*(i+1/2)*x ).*cos( pi*(j+1/2)*y );
            end
        end
        % orthogonalisation step
        [U,R] = qr(A);
    case 'wavelet'
        % wavelet basis
    case 'alpert'
        % alpert basis
        alpert_vm = [2 2]; % piecewise linear
        options.part_type = '2axis'; % dichotomic grouping
        for i=1:n
            v = zeros(n,1); v(i) = 1;  % dirac basis
            [w,info] = perform_alpert_transform_2d(v,xy',alpert_vm, -1, options);
            A(:,i) = w(:);
        end
        % be sure that low frequencies comes first
        U = A(:,end:-1:1);
end


% display
disp('Displaying eigenvectors.');
p = 3;
clf;
for i=1:p^2
    num = 5*(i-1)+1;
    c = U(:,num);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    subplot(p,p,i);
    options.face_vertex_color = c;
    plot_mesh(vertex,face, options);
    str = ['Eigenvector n°', num2str(num+1)];
    title(str);
    shading interp;
end

disp('Press any key to continue.');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot reconstruction
disp('Displaying reconstruction.');

approx_type = 'lin';    % linear
approx_type = 'nonlin'; % non-linear

% reconstruction 
p = 3;
nbr_max_keep = 10; % percent max
clf;
for i=1:p^2
    subplot(p,p,i);

    keep = 1+round(i*nbr_max_keep/p^2); % number of percent of coefficients kept
    % perform the transform
    vertex2 = (U'*vertex')';
    switch lower(approx_type)
        case 'lin'
            % number of kepts coefficients
            m = round(keep/100*n);
            % remove the last coefs
            % vertex2(:,m+1:end) = 0;
            vertex2(:,m+1:end) = 0;
        case 'nonlin'
            % set threshold
            vnorm = sum(vertex2.^2, 1);
            vnorms = sort(vnorm); vnorms = reverse(vnorms);
            thresh = vnorms( round(keep/100*n) );
            % remove small coefs
            mask = vnorm>thresh;
            for k=1:3
                vertex2(k,:) = vertex2(k,:).*mask;
            end
    end
    % reconstruction
    vertex2 = (U*vertex2')';
  
    plot_mesh(vertex2,face);
    str = [num2str(keep), '% of the coefficients'];
    title(str);
end