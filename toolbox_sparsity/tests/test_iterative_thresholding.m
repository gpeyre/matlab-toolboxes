% test for iterative thresholding for sparse coding

path(path, 'toolbox/');

%%%%%%% random setting %%%%%%%%%%
n = 100;
m = n*5;
p = 1; % number of exemplar to sparse code

% dictionary
D = randn(n,m);
D = D./repmat( sqrt(sum(D.^2)), [n,1] );
% exemplars to code
X0 = zeros(m,p);
smax = 5;
for i=1:p
    s = round( 3+rand*(smax-3) );
    sel = randperm(m); sel = sel(1:s);
    X0(sel,i) = 1; % randn(s,1);
end
Y = D*X0;

z = D'*Y;
lambda = std( z(:) )/100;

% conditioning of the dictionary
[a,s,b] = svd(D); 
mu = 1/max(diag(s))^2;

%%% noiseless inversion
options.niter = 1000;
options.lambda_min = 0;
c = D'*Y; c = c(:);
options.lambda_max = mad(z)/0.6745;
options.thresh_type = 'hard';
options.thresh_type = 'soft';
[Xc,Ec] = perform_iterative_thresholding(D,Y(:,1),options);

S = sum(X0>0); slist = min(S(:)):max(S(:));
for s=slist
    I = find(S==s);
    e(s) = norm( Xc(:,I)-X0(:,I), 'fro' ) / (length(I)*s);
end
plot(slist,e(e>0)); axis tight;

%%% comparison of fixed VS decreasing thresholdings
options.niter = 300;
% fixed threshold
options.lambda_min = lambda;
options.lambda_max = lambda;
[Xa,Ea] = perform_iterative_thresholding(D,Y,options);
% decaying threshold
options.lambda_min = lambda;
options.lambda_max = lambda*10;
[Xb,Eb] = perform_iterative_thresholding(D,Y,options);

clf; 
plot( [Ea(:), Eb(:)] );
legend('Fixed', 'Decaying');

