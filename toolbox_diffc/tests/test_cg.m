% test for conjugate gradient

n = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conjugate Gradient
A = randn(n); % random SDP matrix
A = A'*A;
% rand RHS
x0 = randn(n,1);
y = A*x0;

options.is_sdp = 1;
options.niter_max = 200;
options.epsilon = 1e-20;
[x,err,k] = perform_conjugate_gradient(A,y,options);

plot(log2(err)); axis tight;
title('Decreasing of the error.');
xlabel('nb.iter'); ylabel('log(|y-A*x|)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Bi - Conjugate Gradient (works with non SDP)
[U,R] = qr(randn(n));
A = U*diag(rand(n,1))*U';


% rand RHS
x0 = randn(n,1);
y = A*x0;

options.is_sdp = 0;
options.niter_max = 200;
options.epsilon = 1e-20;
[x,err,k] = perform_conjugate_gradient(A,y,options);

plot(log2(err)); axis tight;
title('Decreasing of the error.');
xlabel('nb.iter'); ylabel('log(|y-A*x|)');

x1 = bicg(A,y,1e-10,200);