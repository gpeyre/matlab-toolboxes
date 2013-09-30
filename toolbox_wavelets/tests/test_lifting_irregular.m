% test for lifting scheme on irregular grid

n = 512;

% eta=0 => regular, eta=1 perfect iregular
eta = .99;
x = rescale(rand(n,1),1-eta,1);
x = rescale(cumsum(x));

f = x;

options.x = x;
options.scaling = 1;
Jmin = 4;
fw = perform_wavelet_transform_irregular(f, Jmin, +1, options);
f1 = perform_wavelet_transform_irregular(fw, Jmin, -1, options);
err = norm(f-f1,'fro')/norm(f);
disp(['Error (should be 0):' num2str(err) '.']);

% plot nice wavelet functions
fw = zeros(n,1);
fw( floor(1.5*2^Jmin) ) = 1;
phi = perform_wavelet_transform_irregular(fw, Jmin, -1, options);

clf;
h = plot(x,phi); axis tight;
set(h, 'LineWidth', 2);
