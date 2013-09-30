% test for TV projection on the TV ball

path(path, 'toolbox/');

clear options;
options.bound = 'sym';

rep = 'results/tv-projection/';
if not(exist(rep))
    mkdir(rep);
end
repiter = [rep 'iter/'];
if not(exist(repiter))
    mkdir(repiter);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometrical image
n = 256;
J = 4;
sigma = 100;
options.bound = 'per';
A = ones(n);
for j=0:J-1
    B = A;
    for k=1:2^j
        I = find(B==k);
        U = perform_blurring(randn(n),sigma,options);
        s = median(U(I));
        I1 = find( (B==k) & (U>s) );
        I2 = find( (B==k) & (U<=s) );
        A(I1) = 2*k-1;
        A(I2) = 2*k;
    end
    clf; imageplot(A);
    imwrite(rescale(A), [rep 'geometrical-' num2str(j) '.png'], 'png');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D test
name = 'random';
n = 1024*4;
f = randn(n,1);
tv = compute_total_variation(f, options);

tv_divide = [2 4  8 16 32];
m = length(tv_divide);
clf;
imageplot(f, 'Original', m+1,1,1);
for i=1:m
    tau = tv/tv_divide(i);
    [f1,err_tv,err_l2] = perform_tv_projection(f,tau,options);
    tv1 = compute_total_variation(f1, options);
    imageplot(f1, ['TV/' num2str(tv_divide(i))], m+1,1,i+1);
end
saveas(gcf, [rep name '-1d-tv-projection.png'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load image
name = 'boat';
name = 'noise';
n = 200;
M = rescale( crop( load_image(name), n ) );

sigma = 0; % smoothing
if sigma>0
    M = perform_blurring(M,sigma);
end

% reduce TV norm by projection
tv = compute_total_variation(M, options);
options.niter = 2000;

tv_divide = [2 4 8 16 32 64 128];
m = length(tv_divide);
a = 2; b = ceil((m+1)/2);
options.x = [];
clf;
if strcmp(name, 'noise')
    MM = clamp( 0.45*(M-mean(M(:)))/std(M(:)), -1,1);
else
    MM = M;
end
warning off;
imwrite( rescale(MM), [repiter name '-2d-tv-projection-sig' num2str(sigma) '-01.png'] );
warning on;
imageplot(M, 'Original', a,b,1);
for i=1:m
    tau = tv/tv_divide(i);
    [M1,err_tv,err_l2] = perform_tv_projection(M,tau,options);
    tv1 = compute_total_variation(M1, options);
    disp( ['Final TV error: ' num2str( (tv1-tau)/tv1 ) '.'] );
    options.x = M1;
    if strcmp(name, 'noise')
        MM = clamp( 0.45*(M1-mean(M1(:)))/std(M1(:)), -1,1);
    else
        MM = M1;
    end
    imageplot(MM, ['TV/' num2str(tv_divide(i))], a,b,i+1);
    warning off;
    imwrite( rescale(MM), [repiter name '-2d-tv-projection-sig' num2str(sigma) '-' num2string_fixeddigit(tv_divide(i),2) '.png'] );
    warning on;
end
saveas(gcf, [rep name '-2d-tv-projection-sig' num2str(sigma) '.png'], 'png');
