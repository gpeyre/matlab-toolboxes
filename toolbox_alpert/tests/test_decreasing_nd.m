% test_decreasing_nd - test the nD alpert transform
%
%   Copyright (c) 2004 Gabriel Peyré

alpha = 3;      % regularity of the vertical function.
k = 2;          % number of vanishing moments of the AT.

P = 500;       % number of sampled points
d = 3;         % number of dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position of sampling
pos = rand(d,P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate dD signal
X = 2*(pos - 0.5);
v = cos(2*sum(X.^2)); v = v(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform transform
w = perform_alpert_transform_nd(v,pos,k);

err = l2error(w);
% reglin
fiton = 1:floor(0.2*P);
A = (1:P)';
B = err;
aa = log2(A(fiton));
bb = log2(B(fiton));
[coef,S] = polyfit(aa,bb,1);
disp( sprintf('--> 2D Alpert transform, decreasing=%.4f', coef(1)) );
reglin = log2(A)*coef(1)+coef(2);

% title
thetitle = sprintf('\\alpha=%d, k=%d, decreasing=%.4f', alpha, k, coef(1) );

% plot
clf;
hold on;
plot(log2(1:P), log2(err));
axis tight;
axis([0,log2(P), -10, mmax(log2(err))]);
title(thetitle);

% plot(log2(1:P), log2(errmin), 'k:');
plot(log2(1:P), reglin, 'k:');
hold off;