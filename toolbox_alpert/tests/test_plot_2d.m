% test_plot_2d - test construction of 2D alpert basis
%
%   Copyright (c) 2004 Gabriel Peyré

sampling_type = 'regular';
sampling_type = 'irregular';
% position on a rectangular grid

if strcmp(sampling_type, 'regular')
    k = 4;
    n = 1024;
    
    nx = sqrt(n);
    ny = sqrt(n);
    [Y,X] = meshgrid( 0:1/(nx-1):1, 0:1/(ny-1):1 );
    pos = [X(:)'; Y(:)'];
else 
    k = 4;
    n = 1000;
    
    pos = rand(2,n);
    [x,I] = sort(pos(1,:));
    pos(1,:) = pos(1,I);
    pos(2,:) = pos(2,I);
end


%% plot the vectors
nv = 4;
nh = 4;
clf;
jj = 1;
for i=1:nv
for j=1:nh
    w = dirac(n, floor(n-n/2^(j+jj-1)+i) );
    v = perform_alpert_transform_2d(w, pos, k, -1);
    subplot(nv,nh,j+nh*(i-1));
    % sel = 1:n/2^(4-j);
    % plot(x(sel), v(sel));
    if 0 && strcmp(sampling_type, 'regular')
        M = reshape(v, nx,ny);
        imagesc(-abs(M));
        colormap gray;
        axis image;
    else
        plot_scattered(pos, v);
        shading interp;
        colormap default;
    end
end
end
