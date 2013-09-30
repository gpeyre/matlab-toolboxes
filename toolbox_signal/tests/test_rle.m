% test for run length coding

% some signal with long cluster of 0/1
n = 1024*2;
options.alpha = 0.1;
x = load_signal('regular',n,options); 
x = (x-mean(x))>0;

options.rle_coding_mode = 'shannon';
[tmp,nb_shannon] = perform_rle_coding(x, +1, options);
disp(sprintf('Shannon = %.2f', nb_shannon));

options.rle_coding_mode = 'arithmetic';
[tmp,nb_arith] = perform_rle_coding(x, +1, options);
disp(sprintf('Arithmetic = %.2f', nb_arith));

options.rle_coding_mode = 'arithfixed';
[tmp,nb_arithfixed] = perform_rle_coding(x, +1, options);
disp(sprintf('Arithmetic(laplacian) = %.2f', nb_arithfixed));

[tmp,nb_direct] = perform_arithmetic_coding(x, +1);
disp(sprintf('Direct(entropy) = %.2f', nb_direct));

options.rle_coding_mode = 'nocoding';
[tmp,nb_nocode] = perform_rle_coding(x, +1, options);
disp(sprintf('No code = %.2f', nb_nocode));

% test for bijectivity
options.rle_coding_mode = 'nocoding';
stream = perform_rle_coding(x, +1, options);
xx = perform_rle_coding(stream, -1, options);
disp( sprintf('Error (should be 0): %.2f', norme(x-xx)) );

options.rle_coding_mode = 'arithfixed';
stream = perform_rle_coding(x, +1, options);
xx = perform_rle_coding(stream, -1, options);
disp( sprintf('Error (should be 0): %.2f', norme(x-xx)) );