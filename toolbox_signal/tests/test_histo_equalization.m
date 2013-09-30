% test for equalization of histograms

path(path, 'toolbox/');

name1 = 'lena';
name2 = 'barb';
n = 512;
M1 = rescale( load_image(name1,n) );
M2 = rescale( load_image(name2,n/2) );

M3 = perform_histogram_equalization(M1,M2);
options.nb_bins = 256;
% M3 = perform_histogram_matching(M1,M2);

clf;
imageplot(M1, 'Original', 2,3,1);
imageplot(M2, 'Target', 2,3,2);
imageplot(M3, 'Remapped', 2,3,3); 
subplot(2,3,4);
hist(M1(:),80);
subplot(2,3,5);
hist(M2(:),80);
subplot(2,3,6);
hist(M3(:),80);