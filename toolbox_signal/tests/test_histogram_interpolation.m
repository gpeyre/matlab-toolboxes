% test for histogram interpolation

path(path,'toolbox/');
%M1 = load_image('lena');
%M2 = load_image('barb');

ntest = 2;

n = 50^2;
nbins = 30;

switch ntest
    case 1
        M1 = rescale(randn(n,1),0,1);
        M2 = rand(n,1);
    case 2
        delta = 3;eta = .7;
        M1 = [eta*randn(n/2,1)-delta; eta*randn(n/2,1)+delta];
        M2 = randn(n,1);
    case 3;
        M1 = randn(n,1);
        M2 = randn(n,1)+10;
    case 4;
        M1 = rand(n,1);
        M2 = [rand(n/2,1)*.5-.25; rand(n/2,1)*.5+.75];
    case 5;
        M1 = rand(n,1)-.5;
        M2 = rand(n,1)+rand(n,1)-1;
end

a = min(min(M1(:)),min(M2(:)));
b = max(max(M1(:)),max(M2(:)));

h1 = hist(M1(:),nbins);
h2 = hist(M2(:),nbins);
m = max(max(h1),max(h2));

niter = 100;
for i=1:niter
    options.histinterp = i/niter;
    M = perform_histogram_equalization(M1,M2,options);
    clf;
    hist(M(:),nbins);
    axis([a,b,0,m]);
    drawnow;
end