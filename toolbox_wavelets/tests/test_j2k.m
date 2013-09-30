% test for jpeg2k compression
%
%   Copyright (c) 2005 Gabriel Peyr?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('name')
    name = 'barb';
end

n = 512;
n = [];
M = load_image(name);
M = double( mean(M,3) );
M = crop( M,n);
n = size(M,1);

% options for the wavelets transform
Jmin = 3;
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;

disp('Performing wavelet transform.');
MW = perform_wavelet_transform(M,Jmin,1, options);

bit_depth = 8;
true_comp_ratio = 12;

if ~exist('bpp')
    bpp = 0.5;
end
nbr_bits = (n^2)*bpp;

disp('Performing JPEG2k degradation.');
MW1 = perform_jp2k_degradation(MW,Jmin,nbr_bits,M);

% same but in 2 time
disp('Performing JPEG2k coding/decoding.');
options.Jmin = Jmin;
options.nbr_bits = nbr_bits;
% for i=1:20
stream = perform_wavelet_jp2k_coding(MW, options);
% end
MW2 = perform_wavelet_jp2k_coding(stream, options);
% should be zero
disp( ['Should be 0 : ' num2str(norme(MW1-MW2),2)] );

disp('Performing inverse wavelet transform.');
M1 = perform_wavelet_transform(MW1,Jmin, -1, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display
imageplot( clamp(M/255), 'Original image', 1,2,1 );
imageplot( clamp(M1/255), ['Compressed @' num2str(bpp, 2) 'bpp'], 1,2,2);
