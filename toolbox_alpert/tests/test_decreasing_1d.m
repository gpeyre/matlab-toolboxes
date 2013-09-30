% test_decreasing_1d - test the 2D alpert transform for a C^alpha signal.
% If the Alpert basis has more than alpha vanishing moments
% (check the value of k), then the decreasing of the error
% of reconstruction with M coefficients should be like
% |f-f_M| ~ M^{-alpha}
%
%   Copyright (c) 2004 Gabriel Peyré

alpha = 3;      % regularity of the function.
k = 3;          % number of vanishing moments of the AT.

% generate a super resolution signal
nn = 5000;

a = gen_signal(nn, alpha);

% number of points
P = 600;

% sample at random location
x = floor( rand(1,P)*nn ) + 1;
v = a(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform transform
tic;
clear options;
options.use_mex= 1;
w = perform_alpert_transform_1d(v,x,k, 1, options);
disp(['1d using mex takes ' num2str(toc)]);
%%%%
tic;
options.use_mex= 0;
w1 = perform_alpert_transform_1d(v,x,k, 1, options);
disp(['1d not using mex takes ' num2str(toc)]);
%%%%
tic;
V = build_alpert_matrix_1d(x,k);
w2 = (V')*v;
disp(['1d matrix takes ' num2str(toc)]);
disp(['difference between the two (should be 0) : ' num2str(norme(w2-w))]);    % should be 0

err = l2error(w);
err1 = l2error(w1);

% reglin
fiton = 1:floor(P/4);
A = (1:P)';
B = err;
aa = log2(A(fiton));
bb = log2(B(fiton));
[coef,S] = polyfit(aa,bb,1);
disp( sprintf('----> 1D Alpert transform, decreasing=%.4f', coef(1)) );
reglin = log2(A)*coef(1)+coef(2);


% plot
plot_on = 1:floor(P/2);
clf;
plot(log2(plot_on), log2(err(plot_on)), 'b', log2(plot_on), log2(err1(plot_on)), 'r.', log2(plot_on), reglin(plot_on), 'b:');
legend('using mex', 'not using mex', 'reglin');
axis tight;
% axis([0,log2(P), -10, max(max(log2(err)))]);
title(sprintf('1D alpert transform, \\alpha=%d, decreasing=%.4f', alpha, coef(1)));