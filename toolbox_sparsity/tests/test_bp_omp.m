% test omp vs bp

n = 100;
m = 200;
s = 20;
p = 100; % nbr tests

err_bp = [];
err_omp = [];

D = randn(n,m);
D = D./repmat( sqrt(sum(D.^2)), [n 1] );

X0 = randn(m,p);
t = sort(abs(X0));
t = t(end-s,:);
t = repmat(t, [m 1]);
X0 = X0 .* (abs(X0)>t);
Y = D*X0;

options.use_mex = 1;
for k=1:2*s
    progressbar(k,2*s);
    options.nbr_max_atoms = k;
    options.sparse_coder = 'mp';
    X = perform_omp(D,Y,options);
    err_bp(end+1) = norm(Y-D*X,'fro');
    
    options.sparse_coder = 'omp';
    X = perform_omp(D,Y,options);
    err_omp(end+1) = norm(Y-D*X,'fro');
end

clf;
plot(1:2*s, err_omp, 1:2*s, err_bp);
axis tight;
legend('OMP', 'BP');
