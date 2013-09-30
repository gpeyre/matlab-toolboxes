% test for the quincunx wavelet transform

name = 'barb';
name = 'lena';
n = 512;
M = load_image(name, n);
M = rescale(M);

Jmax = log2(n)-1;
Jmin = Jmax-2;

rep = 'results/quincunx/';
if not(exist(rep))
    mkdir(rep);
end

options.null = 0;
% Forward transform
fprintf('Computing forward transform ... ');
[MW,options.quincunx_filters] = perform_quincunx_wavelet_transform(M,Jmin,+1,options);
fprintf('done.\n');
% Backward transform
fprintf('Computing backward transform ... ');
M1 = perform_quincunx_wavelet_transform(MW,Jmin,-1,options);
fprintf('done.\n');


% error of reconstruction
disp([ 'Error of reconstruction ' num2str(psnr(M,M1), 4) 'dB.']);

% display
options.gamma = .75;
clf;
MW1 = plot_quincunx_wavelet(MW, Jmin, options);
warning off;
imwrite(MW1(1:end/2,1:end/2), [rep name '-quincunx-wav.png'], 'png');
imwrite(rescale(M), [rep name '-original.png'], 'png');
warning on;


% select some sub-images
options.transform = 'quincunx';
m = 2*(Jmax-Jmin+1)+1;
nb = ceil(m/2);
bnd = 1;
k = 0;
clf;
for j=Jmax:-1:Jmin
    qq = 1:2;
    if j==Jmin
        qq(end+1) = 0;
    end
    for q=qq
        k = k+1;
        [selx,sely] = compute_quadrant_selection(j,q, options);
        MWj = MW(selx,sely);
        
        MWj = MWj/max(abs(MWj(:)));
        % gamma boosting of contrast
        if q==0
            MWj = rescale(MWj);
        else
            MWj = abs(MWj.^options.gamma).*sign(MWj);
            MWj = (MWj+1)/2;
            MWj(1,:) = bnd; MWj(end,:) = bnd; MWj(:,1) = bnd; MWj(:,end) = bnd;
        end
        
        if q==1
            nj = size(MWj,2);
            MH = zeros(2*nj,2*nj)+1;
            MH(1:2:end,1:2:end) = MWj(1:2:end,:);
            MH(2:2:end,2:2:end) = MWj(2:2:end,:);
            [Y,X] = meshgrid(1:2*nj,1:2*nj);
            % rotate
            U = X+Y; V = -X+Y+nj*2;
            T = ones(4*nj);
            T(U+(V-1)*4*nj) = MH;
            T = T(2:2:end,2:2:end);
        else
            T = MWj;
        end
        
        subplot(2,nb,k);
        imagesc(T);
        title(['j=' num2str(j) ', q=' num2str(q)]);
        axis image; axis off;
        
        warning off;
        imwrite(T, [rep name '-quincunx-wav-j' num2str(j) '-q' num2str(q) '.png'], 'png');
        warning on;
    end
end
colormap gray(256);