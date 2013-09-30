% test for low rank approximation of an image

name = 'lena';
n = 128;
M = load_image(name);
M = crop(M,n);

[U,S,V] = svd(M);


%% display reconstruction
m_list = [5 10 50 100];
for i = 1:length(m_list)
    m = m_list(i);
    s = diag(S);
    s(m+1:end) = 0;
    M1 = U*diag(s)*V';
    subplot(2,2,i);
    imagesc(M1); axis image; axis off;
    title([num2str(m) ' eigv']);
end
colormap gray(256);

return;

% test for pseudo inverse
n = 300; m = 4*n;
D = rand(n,m);
tic;
B = D' * (D*D')^(-1);
t = toc;
tic;
B1 = pinv(D);
t1 = toc;
e = norm(B-B1, 'fro')/norm(B, 'fro');
disp(['speed-up=' num2str((t1-t)/t*100) '%, err=' num2str(e*100)]);