% test for homotopy solving

path(path, 'toolbox/');

% number of measurements
n = 100;
% number of atoms
m = 1000;
% sparsity
s = 15;

% random signal
x = compute_rand_sparse(m,s);
% random measurements matrix
P = compute_compressed_sensing_matrix(n,m);
% measures
y = P*x;

%% charles homotopy
gamma = 1e-9;
Xhomot = Homotopie2(y,P,gamma,0,x);

%% candes basis pursuit
x0 = x; % x0 = zeros(m,1);
Xcandes = l1qc_logbarrier(x,P,P',y,0.01);

%% sparselab BP
lambda = 0; maxIters = 100; OptTol = 1e-5;
tic;
Xspars = SolveBP(P, y, m, maxIters, lambda, OptTol);
disp(['Sparselab timing: ' num2str(toc) 'sec']);


%% sparselab Lasso
algType = 'lars';
verbose = 0;
lambdaStop = 0;
solFreq = 0;
maxIters = 1000;
Xlasso = SolveLasso(P, y, m, algType, maxIters, lambdaStop, solFreq, verbose, OptTol);

%% sparselab OMP
maxIters = 1000;
Xomp = SolveOMP(P, y, m, maxIters, lambdaStop, solFreq, verbose, OptTol);

%% Boyd BP
%run the l1-regularized least squares solver
tic;
lambda = 0.01;
rel_tol = 0.01;
[Xboyd,status]=l1_ls(P,P',n,m,y,lambda,rel_tol);
disp(['Boyd BP timing: ' num2str(toc) 'sec']);


Fcoef = CalculF(P,x);
ERCcoef = 1-ercm2(P,x);
disp( ['Fuchs=' num2str(Fcoef) ', ERC=' num2str(ERCcoef)] );

sol = {Xhomot, Xcandes, Xspars, Xlasso, Xomp};
lgd = {'Charles homot', 'BP (Candes)', 'BP (Sparselab)', 'Lasso', 'OMP'};

p = ceil(length(sol)/2);
clf;
disp(['--> Original L1=' num2str( sum( abs(x) ) ) ]);
ax = [];
for i=1:length(sol)
    x1 = sol{i};
    err = norm(x1-x, 'fro');
    L1 = sum( abs(x1) );
    str = [lgd{i} ', L1= ' num2str(L1) ', err=' num2str(err)];
    disp(['--> ' str]);
    ax(end+1) = subplot(2,p,i);
    stem(x1); axis tight; 
    title(str);
end
linkaxes(ax, 'xy');

