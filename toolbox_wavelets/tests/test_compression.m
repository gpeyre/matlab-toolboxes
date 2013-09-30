%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for wavelet based compression.
%
% This script test various kind of wavelet coefficient encoding.
% A compressor is based on the following :
%   * a transform (here a 9/7 wavelet biorthogonal transform).
%   * a quantizer (here either a uniform quantizer for arithmetic 
%		or a bit-plane encoding for SPIHT).
%	* a binary coding scheme (either a simple arithmetic encoding
%		of the set of quantized coefficient or a SPIHT encoding).
%
%   Copyright (c) 2006 Gabriel Peyré
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path, 'toolbox/');

name = 'barb';
name = 'lena';

M = load_image(name);

% reduce size for speed up
n = 256;
M = M(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2);
M = rescale(M,0,255);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a set of threshold (valid for a 0...255 valued image)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a set of threshold (valid for a 0...255 valued image)
%  T=10 .. 50   ==>  bpp=1 .. 0.2
%  T=2.2 .. 10  ==>  bpp=3 .. 1
%  T=0.5 .. 2.2 ==>  bpp=5.5 .. 3
T_high = [8:2:20 20:5:50]; % for high compression ranges 50:10:80
T_low  = [2:0.6:5 6:1.2:10]; % 0.5:0.06:0.5 0.5:0.3:2 
T_high = [8:3:20 20:8:50]; % for high compression ranges 50:10:80
T_high = [5:5:20 20:12:40]; % for high compression ranges 50:10:80
T_low  = [2:0.6:5 6:1.3:10]; % 0.5:0.06:0.5 0.5:0.3:2 
T_range = 'high';
if strcmp(T_range, 'low')
    T = T_low;
elseif strcmp(T_range, 'high')
    T = T_high;
else
    T = [T_low T_high];
end
T = sort( unique(T) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;
Jmin = 3;
options.Jmin = Jmin;
MW = perform_wavelet_transform(M,Jmin, 1, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to reccord the results
psnr_spiht = [];
bit_spiht = [];
psnr_jp2k = [];
bit_jp2k = [];

test_spiht = 0;
test_jp2k = 1;
tested_artihmetic = [5 1 1 2 2 7];
tested_mode =       [1 1 2 1 2 1];
lgd_arithmetic = {'Shannon', 'LetItWave', 'LetItWave packed', 'Escape', 'Escape Packed', 'Gen.Laplacian'};
for k=1:length(tested_artihmetic)
    psnr_arith{k} = [];
    bit_arith{k} = [];
end

% test for psnr
fprintf( 'Testing for T=' );
for i = 1:length(T)
    t = T(i);
    % direct transform
    fprintf( sprintf('%.1f,', t) );
    
    % arithmetic coding
    for k=1:length(tested_artihmetic)
        options.coder_type = tested_artihmetic(k);
        options.coding_mode = tested_mode(k);
        [stream,R_arith] = perform_wavelet_arithmetic_coding(MW, t, options);  % coding
        MW_arith = perform_wavelet_arithmetic_coding(stream, t, options);      % decoding
        M_arith = perform_wavelet_transform(MW_arith, Jmin, -1, options);   % untransform
        psnr_arith{k} = [psnr_arith{k}, psnr( M_arith,M )];
        bit_arith{k} = [bit_arith{k}, R_arith/n^2];
    end
    
    % spiht coding
    if test_spiht
        options.Jmin = Jmin;      % minimum scale of the transform
        options.nb_bits = R_arith;   % target number of bits
        [stream,R_spiht] = perform_spiht_coding(MW,options);                % coding
        MW_spiht = perform_spiht_coding(stream);                            % decoding
        M_spiht = perform_wavelet_transform(MW_spiht, Jmin, -1, options);   % untransform
        psnr_spiht = [psnr_spiht, psnr( M_spiht,M )];    
        bit_spiht = [bit_spiht, R_spiht/n^2];
    end
    
    % jpg2k coding
    if test_jp2k
        [MW_jp2k,R_jp2k] = perform_jp2k_degradation(MW,Jmin,R_arith,M);
        M_jp2k = perform_wavelet_transform(MW_jp2k, Jmin, -1, options);   % untransform
        psnr_jp2k = [psnr_jp2k, psnr( M_jp2k,M )];    
        bit_jp2k = [bit_jp2k, R_jp2k/n^2];
    end
end
fprintf('\n');

x = []; y = []; lgd = {};
for i=1:length(tested_artihmetic)
    lgd{end+1} = ['Arithmetic ' lgd_arithmetic{i}];
    x = [x bit_arith{i}(:)];
    y = [y psnr_arith{i}(:)];
end
if test_spiht
    lgd{end+1} = 'SPIHT';
    x = [x bit_spiht(:)];
    y = [y psnr_spiht(:)];
end
if test_jp2k
    lgd{end+1} = 'JPEG2k';
    x = [x bit_jp2k(:)];
    y = [y psnr_jp2k(:)];
end
plot(x, y);
axis tight;
legend(lgd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute some nice plot for psnr gain 
m1 = min( x(:) );
m2 = max( x(:) );
bpp = linspace(m1,m2,200);
bpp = sort( unique([bpp(:); x(:)]) );

interp_psnr = {};
for i=1:size(x,2)
    interp_psnr{i} = interp1(x(:,i), y(:,i), bpp);
end
gain_x = {};
gain_y = {};
for i=2:size(x,2)
    gain_y{i-1} = interp_psnr{i}-interp_psnr{1};
    I = find(~isnan(gain_y{i-1}));
    gain_y{i-1} = gain_y{i-1}(I);
    gain_x{i-1} = bpp(I);
end

clf;
hold on;
for i=1:length(gain_x)
    plot(gain_x{i}, gain_y{i}, get_color_from_index(i));
end
title(['PSNR Gain with respect to ' lgd{1}]);
xlabel('bbp');
ylabel('+PSNR');
axis tight;
legend({lgd{2:end}});
hold off;