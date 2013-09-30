% test of Laplacian pyramid

n = 256;
name = 'lena';
M = load_image(name,n);
Jmin = 3;
MW = perform_pyramid_transform_simoncelli(M, Jmin);

save_image = 0;

for i=1:4
    subplot(2,2,i);
    imagesc(MW{i});
    axis image;
    axis off;
    if save_image
        imwrite(rescale(MW{i}), sprintf('MW%d.png', i), 'png');
    end
end
colormap gray(256);

if save_image
    imwrite(rescale(M), sprintf('M.png', i), 'png');
end