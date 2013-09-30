% test for texture synthesis

path(path, 'images/');
path(path, '../images/');

name = 'grass';
name = 'turbulence1';
name = 'fingerprint3';
name = 'dunes';

n = 128;
M = load_image(name, n);
M = crop(M,n);

s = size(M,3);
n1 = 128;
M1 = randn(n1,n1,s);

options.niter_synthesis = 5;
options.verb = 1;
disp('--> Steerable synthesis.');
options.synthesis_method = 'steerable';
Msteer = perform_wavelet_matching(M1,M,options);
disp('--> Wavelet orthogonal synthesis.');
options.synthesis_method = 'wavelets-ortho';
Morth = perform_wavelet_matching(M1,M,options);
disp('--> Quincunx synthesis.');
options.synthesis_method = 'quincunx-ti';
Mquin = perform_wavelet_matching(M1,M,options);
disp('--> Wavelets TI synthesis.');
options.synthesis_method = 'wavelets-ti';
Mwavti = perform_wavelet_matching(M1,M,options);

clf;
imageplot({M Msteer Morth Mquin Mwavti},{'Original', 'Synth. steerable', 'Synth. wav. ortho', 'Synth. quincunx', 'Synth. wav. ti'});