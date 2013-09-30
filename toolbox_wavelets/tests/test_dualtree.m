% test for dual tree transform

name = 'lena';
name = 'disk';

n = 256;
M = load_image(name,n);
Jmin = 4;
MWr = perform_real_dualtree_transform(M,Jmin);
MWc = perform_cpx_dualtree_transform(M,Jmin);


clf;
for i=1:6
    subplot(2,3,i);
    imagesc(MWr{i});
    axis image; axis off;
end
colormap gray(256);