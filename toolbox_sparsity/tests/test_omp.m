% test for orthogonal matching pursuit procedure
n = 200; % dimensionality of the problem
p = 1000; % number of basis elements

% random dictionary
D = randn(n,p);
D = D ./ repmat( sqrt(sum(D.^2,1)), [n 1] );

% generate sparse vector
k = 30;
s = 50; % number of tests
x0 = zeros(p,s);
for i=1:s
    pos = randperm(p); pos = pos(1:k);
    x0( pos,i ) = 1;
end
f = D*x0;

options.tol = 1e-3;  % tolerance
options.nbr_max_atoms = 200;
options.nbr_max_atoms = k;

% perform OMP
if 0
tic;
options.use_slow_code = 1;
x = perform_omp(D,f,options);
toc;
end

tic;
options.use_slow_code = 0;
x1 = perform_omp(D,f,options);
disp( ['Matlab code: ' num2str(toc) 'sec.'] );

tic;
options.use_mex=1;
x2 = perform_omp(D,f,options);
options.use_mex=0;
disp( ['Mex code:    ' num2str(toc) 'sec.'] );


%% display L1 norm and error
clf; subplot(3,1,1);
plot(1:s, sum(abs(x1),1), 1:s, sum(abs(x2),1));
axis tight;
xlabel('Trial #'); ylabel('L1 norm');
legend('Matlab', 'Mex');

e1 = sqrt( sum( (D*x1-f).^2, 1 )/n );
e2 = sqrt( sum( (D*x2-f).^2, 1 )/n );

subplot(3,1,2);
plot(1:s, e1, 1:s, e2);
axis tight;
xlabel('Trial #'); ylabel('Error |A*x-A*x0|^2');
legend('Matlab', 'Mex');

e1 = sqrt( sum( (x1-x0).^2, 1 )/k );
e2 = sqrt( sum( (x2-x0).^2, 1 )/k );

subplot(3,1,3);
plot(1:s, e1, 1:s, e2);
axis tight;
xlabel('Trial #'); ylabel('Error |x-x0|');
legend('Matlab', 'Mex');

