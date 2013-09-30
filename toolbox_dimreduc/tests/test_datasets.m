% test for loading of face dataset

name_list = {'umist_cropped','frey_rawface','olivettifaces'};
options.nclass=1:20;


for i=1:length(name_list)
    p = [20 20];
    name = name_list{i};
    options.nbr = prod(p);
    M = load_images_dataset(name, options);
    nbr = size(M,3);
    p = min( p, ceil([sqrt(nbr) sqrt(nbr)]) );
    s=0;
    n = [size(M,1) size(M,2)];
    p = min( p, round([2000 2000]./n) );
    A = zeros(p.*n);
    for kx=1:p(1)
        for ky=1:p(2)
            s = s+1;
            if s<nbr
                selx = (kx-1)*n(1)+1:kx*n(1);
                sely = (ky-1)*n(2)+1:ky*n(2);
                A(selx,sely) = M(:,:,s);
            end
        end
    end
    figure;
    imagesc(A);
    axis image; axis off;
    colormap gray(256);
end