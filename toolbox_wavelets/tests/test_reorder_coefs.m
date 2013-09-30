% test for coefficient reordering Mallat<->Inplace

n = 256;
% only Jmin==0 is accepted for now
Jmin = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D test
x = 1:n;
y = perform_haar_transform(x,Jmin, 1);
y1 = reorder_coefs(y, Jmin,1);

subplot(2,2,1);
plot(x); axis tight;
title('Original signal.');

subplot(2,2,2);
plot(y); axis tight;
title('Inplace transformed.');

subplot(2,2,3);
plot(y1); axis tight;
title('Mallat''s order.');

yy = reorder_coefs(y1, Jmin ,-1);
disp( sprintf('Reordering Mallat->inplace error (should be 0): %.2f', norm(yy-y, 'fro')) );
yy = reorder_coefs(y, Jmin ,1);
disp( sprintf('Reordering inplace->Mallat error (should be 0): %.2f', norm(yy-y1, 'fro')) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D test
x = load_image('lena', n);
y = perform_haar_transform(x,Jmin, 1);
y1 = reorder_coefs(y, Jmin, 1);

clf;
subplot(1,2,1);
imagesc(x); axis image; axis off;
title('Original signal.');

subplot(1,2,2);
plot_wavelet(y1,Jmin); axis image;
title('Mallat''s order.');

yy = reorder_coefs(y1, Jmin ,-1);
disp( sprintf('Reordering Mallat->inplace error (should be 0): %.2f', norm(yy-y, 'fro')) );
yy = reorder_coefs(y, Jmin ,1);
disp( sprintf('Reordering inplace->Mallat error (should be 0): %.2f', norm(yy-y1, 'fro')) );