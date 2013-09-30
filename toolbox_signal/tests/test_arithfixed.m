% test for fixed length arithmetic coding

if ~exist('x')
    n = 2048;
    x = rand(n,1).^4;
    q = 4; % number of tokens
    x = floor(x*q)+1;
end

options.histo = compute_histogram(x);

%%% Fast mex code %%%
options.coder_type = 8;
tic;
stream = perform_arithmetic_coding(x, +1, options);
x1 = perform_arithmetic_coding(stream, -1, options);
disp(['Mex Code timing:     ' num2str(toc,2) '.']);
disp(['Error (should be 0): ' num2str(norme(x-x1),2) '.']);
nbr_bits_mex = length(stream)*8;

%%% Slow matlab code %%%
options.coder_type = 6;
tic;
stream1 = perform_arithmetic_coding(x, +1, options);
x1 = perform_arithmetic_coding(stream1, -1, options);
disp(['Matlab Code timing:  ' num2str(toc,2) '.']);
disp(['Error (should be 0): ' num2str(norme(x-x1),2) '.']);
nbr_bits_mat = length(stream1);

nbr_bits_entropy = ceil( compute_entropy(x)*length(x) );

disp(['Mex #bits:     ' num2str(nbr_bits_mex)]);
disp(['Matlab #bits:  ' num2str(nbr_bits_mat)]);
disp(['Entropy #bits: ' num2str(nbr_bits_entropy)]);