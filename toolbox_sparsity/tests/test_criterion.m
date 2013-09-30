% test for criterions

path(path, 'toolbox/');

n = 100;
m = n*10;
s = 16;

D = compute_compressed_sensing_matrix(n,m);

x = compute_rand_sparse(m,s);
y = D*x;

tic;
[f,i,si] = compute_fuchs_criterion(D,x);
toc;
tic;
[f1,i1,si1]=CalculF2(D,x);
toc;

disp(['Should be 0: ' num2str(i-i1), ', ' num2str(si-si1)]);

tic;
niter = 100;
[x1,succeded,flist1] = perform_extension_support(D,x,niter,0);
f1 = min(flist1); % F*(x0)
toc;
tic;
[x2,succeded,flist2] = perform_extension_support_adv(D,x,niter);
f2 = min(flist2);
toc;
[x3,f3] = perform_extension_cleaning(D,x0,x1);
[x4,f4] = perform_extension_cleaning(D,x0,x2);

% solve BP
maxIters = 20; lambda = 0; OptTol = 1e-5;
x_bp = SolveBP(D, y, m, maxIters, lambda, OptTol);
rec_bp = compute_recovery_rate(x,x_bp);


% solve OMP
maxItersOMP = s; solFreq = 0; verbose = 0; lambdaStop = 0;
x_omp = SolveOMP(D, y, m, maxIters, lambdaStop, solFreq, verbose, OptTol);
rec_bp = compute_recovery_rate(x,x_omp);