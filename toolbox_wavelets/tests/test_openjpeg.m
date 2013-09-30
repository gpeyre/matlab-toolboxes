% test for OpengJPEG implementation of the JPEG2000 Codec
clear options;
name = 'lena';
M = load_image(name);

test_type = 'rate';
test_type = 'db';
test_type = 'bits';

switch test_type
    case 'rate'
        % compression of a factor 30
        options.rate = 30;
    case 'db'
        % compression with given erro
        options.db = 30;
    case 'bits'
        options.nbr_bits = 5e4;
end

% Compression:
[stream,nbr_bits] = perform_jp2k_compression(M,options);
% De-compression:
M1 = perform_jp2k_compression(stream,options);
bd = ceil(log2(max(M(:))));
nbr_bits_orig = bd*prod(size(M));

clf;
subplot(1,2,1);
imagesc(M); axis image, axis off;
title('Original');
subplot(1,2,2);
imagesc(M1); axis image, axis off;
title(['Compressed, db=' num2str(psnr(M,M1,2^bd),3) ', rate=' num2str(nbr_bits_orig/nbr_bits,3), ' nbrbits=' num2str(nbr_bits,3) ]);
colormap gray(256)