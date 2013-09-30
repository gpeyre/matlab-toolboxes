% test for the tensor structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path, 'toolbox/');

n = 128;
rep = 'faces/';
rep = '';
name = 'polygons_blurred';
name = 'rakan';
name = 'barb';
name = 'lena';
M0 = load_image([rep name]);
M = sum(M0,3);
M = M(end/2-n/2+1:end/2+n/2,end/2-n/2+1:end/2+n/2);
M0 = M0(end/2-n/2+1:end/2+n/2,end/2-n/2+1:end/2+n/2,:);
M = M + 0.1*randn(n);


sigma1 = 1.5; 
sigma2 = 4;


options.m_theta = 16;

options.use_anisotropic=0;
options.use_renormalization = 0;
disp('Computing isotropic tensors.');
H = compute_structure_tensor(M,sigma1,sigma2, options);

options.use_anisotropic=0;
options.use_renormalization = 1;
options.sigmat = 0.4;
disp('Computing anisotropic tensors.');
H1 = compute_structure_tensor(M,sigma1,sigma2, options);

clf;
subplot(1,2,1);
plot_tensor_field(H, M0, 5);
title('Isotropic');

subplot(1,2,2);
plot_tensor_field(H1, M0, 5);
title('Anisotropic');

rep = 'results/';
if ~exist(rep)
    mkdir(rep);
end
saveas( gcf, [rep name '_tensor.png'], 'png' );

return;

[e1,e2,l1,l2] = perform_tensor_decomp(H1);

an = l1./l2; ae = l1.*l2;
ae_target = 4; ae = ae_target * ae./mean(ae(:));
an_median = 4; 
an = 2.^( log2(an_median)/log2(median(an(:))) * log2(an)   );
l1 = sqrt(ae.*an);
l2 = sqrt(ae./an);
H2 = perform_tensor_recomp(e1,e2,l1,l2);

close all;
plot_tensor_field(H, M0, 5);
figure;
plot_tensor_field(H2, M0, 5);
