% test for ti quicunx

name = 'turbulence';
name = 'lena';
n = 256;
M = load_image(name);
M = rescale( crop(M,n) );
Jmax = log2(n)-1;
Jmin = Jmax-5;

% boundary conditions
options.bound = 'per';
options.bound = 'sym';
% vanishing moments
vm = 6;
options.primal_vm = vm;
options.dual_vm = vm;

% transform
disp('Computing forward transform');
MW = perform_quicunx_wavelet_transform_ti(M,Jmin,options);
disp('Computing backward transform');
M1 = perform_quicunx_wavelet_transform_ti(MW,Jmin,options);
disp(['Error=' num2str(norm(M-M1,'fro')/norm(M,'fro')) ' (should be 0).']);

rep = 'results/quincunx/';
if not(exist(rep))
    mkdir(rep);
end

% display dual wavelets
clf;
for j=7:12
    MW = MW*0;
    MW(end/2,end/2,j) = 1;
    M1 = perform_quicunx_wavelet_transform_ti(MW,Jmin,options);
    imageplot(M1, '', 4,3,j);
    % save
    if j>6
        warning off;
%        imwrite(rescale(M1), [rep 'quincunx-wavelet-' num2string_fixeddigit(j,2) '.png'], 'png');
        warning on;
        clf;
        surf(M1);
        shading interp; colormap jet(256)
        view(-20,50); axis tight; axis off;
        camlight;
        saveas(gcf, [rep 'quincunx-wavelet-vm' num2str(vm) '-' num2string_fixeddigit(j,2) '.png'], 'png');
    end
end