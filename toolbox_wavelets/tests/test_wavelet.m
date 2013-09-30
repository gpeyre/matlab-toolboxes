% Test file for the wavelet library.
%   
%   Copyright (c) 2003 Gabriel Peyré

disp('---> ND Wavelet test.');
M = gen_signal_2d(256,3);
MM = fwt_isotropic(M,2);
MM = fwt_hyperbolic(M,2);

disp('---> 1D wavelet test.');
timing1 = [];
for n=2.^(1:20)
    J = floor( log2(n) );
    f = gen_signal(n,2);
    tic;
    g = perform_79_transform(f, J);
    t = toc; timing1 = [timing1, t];
    disp( sprintf('%d coefficients : %.2f sec.', n, t) );
end