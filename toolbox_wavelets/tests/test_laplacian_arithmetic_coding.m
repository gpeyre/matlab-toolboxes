% test for laplacian arithemtic coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 512;
Jmin = 4;
name = 'lena';
M = load_image(name,n);
MW = perform_wavelet_transform(M,Jmin,+1);


A = MW(1:end/2, end/2+1:end);
n = size(A,1);
% quantize to finite precision
T = 10;
[tmp,A] = perform_quantization(A,T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arithmetic cosing using generalized laplacian model
options.coder_type = 7;
fprintf('Performing coding ... ');
[stream,nbr_bits_lapl] = perform_arithmetic_coding(A(:),+1,options);
fprintf(' done.\n');
fprintf('Performing decoding ...');
A1 = perform_arithmetic_coding(stream,-1,options);
fprintf(' done.\n');
A1 = reshape(A1,n,n);
disp( ['Error (should be 0) ' num2str(norme(A-A1))] );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arithmetic cosing using adaptive coding
options.coder_type = 1;
fprintf('Performing coding ... ');
[stream,nbr_bits_adapt] = perform_arithmetic_coding(A(:),+1,options);
fprintf(' done.\n');
fprintf('Performing decoding ...');
A1 = perform_arithmetic_coding(stream,-1,options);
fprintf(' done.\n');
A1 = reshape(A1,n,n);
disp( ['Error (should be 0) ' num2str(norme(A-A1))] );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shannon bound
nbr_bits_shannon = compute_entropy(A(:))*length(A(:));

disp(['Adaptive:  ' num2str(round(nbr_bits_adapt)) 'bits.']);
disp(['Laplacian: ' num2str(round(nbr_bits_lapl)) 'bits.']);
disp(['Shannon:   ' num2str(round(nbr_bits_shannon)) 'bits.']);