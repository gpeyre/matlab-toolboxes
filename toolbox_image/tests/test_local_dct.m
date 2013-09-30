% test for the DCT2 and DCT4 transform

name = 'lena';

n = 256;
M = load_image(name);
M = crop(M,n);

M = randn(n);

w = 16;
q = 8;

options.window_type = 'sin';

%% orthogonal transforms
% local overlapping orthogonal DCT4 transform
tic;
MD1 = perform_windowed_dct4_transform(M, w, +1, options);
disp(['Time for Local DCT4: ' num2str(toc) 'sec.']);
% local non-overlapping orthogonal DCT2 transform
tic;
MD2 = perform_local_dct_transform(M, +1, w);
disp(['Time for orthogonal DCT2: ' num2str(toc) 'sec.']);



%% redundant transforms
% local overlapping redundant transform
tic;
MD3 = perform_windowed_dct_transform(M, w, q, n, options);
disp(['Time for Local DCT2: ' num2str(toc) 'sec.']);
% local overlapping redundant fourier transform (aka Gabor)
tic;
MF = perform_windowed_fourier_transform(M, w,q,n, options);
disp(['Time for Local FFT: ' num2str(toc) 'sec.']);

% check energy conservation, should be tight frame expansion
e = sqrt( sum(M(:).^2)/n^2 );
ef = sqrt( sum(abs(MF(:)).^2)/n^2 );
ed1 = sqrt( sum(abs(MD1(:)).^2)/n^2 );
ed2 = sqrt( sum(abs(MD2(:)).^2)/n^2 );
ed3 = sqrt( sum(abs(MD3(:)).^2)/n^2 );
disp(['Energies (should be the same): ' num2str(e) ', ' num2str(ef) ', ' num2str(ed1) ', ' num2str(ed2) ', ' num2str(ed3)  '.' ]);

% test for alternate normalization with unit norm basis vectors
options.normalization = 'unit';
tic;
MD4 = perform_windowed_dct_transform(M, w,q,n, options);
disp(['Time for Local DCT2+renorm: ' num2str(toc) 'sec.']);
MM = perform_windowed_dct_transform(MD4, w,q,n, options);
options.normalization = 'tightframe';
disp(['Error (should be 0): ' num2str(norm(M-MM)) '.']);
ed4 = sqrt( sum(abs(MD4(:)).^2)/n^2 );

% check to see if basis vectors have unit norm, since we start from white noise
std(MD3(:))
std(MD4(:))
