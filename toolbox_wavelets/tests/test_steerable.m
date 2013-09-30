% test for steerable transform

n = 512;
name = 'turbulence';
name = 'barb';
name = 'disk';
name = 'lena';
M = load_image(name,n);

k = 4;
options.nb_orientations = k;
J = 4;
Jmax = log2(n)-1;
Jmin = Jmax-J+1;
MW = perform_steerable_transform(M,Jmin,options);
M2 = perform_steerable_transform(MW,Jmin,options);

save_image = 1;
rep = ['results/steerable/'];

% reconstruction error
disp(['--> Reconstruction error: ' num2str(psnr(M,M2)) 'dB']);

if ~exist(rep)
    mkdir(rep);
end


% display and save the images
m = 0;
clf;
warning off;
for j=1:J
    for s=1:k
        m = m+1;
        A = MW{m+1};
        % A = rescale_wavelet_coefs(M,eta);
        str = ['j' num2str(j) '-s' num2str(s)];
        imageplot(A, str, k,J,m)
        if save_image
            imwrite(rescale(MW{m+1}), [rep name '-' str '.png'], 'png');        
%            imwrite(rescale(MW{m+1}), [rep name '-' str '.jpg'], 'jpg');   
        end
    end
end
colormap gray(256);


if save_image
    imwrite(rescale(M), [rep name '-original.png'], 'png');
    imwrite(rescale(MW{1}), [rep name '-high.png'], 'png');
    imwrite(rescale(MW{end}), [rep name '-low.png'], 'png');
end
warning on;