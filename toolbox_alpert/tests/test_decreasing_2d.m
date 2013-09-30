% test_decreasing_2d - test the 2D alpert transform
%
% If the Alpert basis has more than alpha vanishing moments
% (check the value of k), then the decreasing of the error
% of reconstruction with M coefficients should be like
% |f-f_M| ~ M^{-alpha/2}
%
% You have to remember that this is not a fully
% 2D wavelet construction (the transform
% is pyramidal only on 1 dimension, on 
% the other this is just a fixed polynomial basis
% with a given number of slices).
%
% Here we plot the decreasing of the error for various number of 
% slices. You can see that as the number of coefficient M 
% increases, the number of slices increases.
% This is because you need more precision also on 
% the Y direction to reconstruct the function
% (which has the same regularity both in the X and in the Y direction).
%
%   Copyright (c) 2004 Gabriel Peyré

clear;
alpha = 3;      % regularity of the vertical function.
k = 4;          % number of vanishing moments of the AT.

P = 30^2;       % number of sampled points


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position of sampling
type_sampling = 'random';

if strcmp(type_sampling,'random')
    % random sampling
    x = rand(1,P);
%    x = sort(x);
    y = rand(1,P);

elseif strcmp(type_sampling,'uniform')
    % uniform distribution
    xx = 0:1/(sqrt(P)-1):1;
    [X,Y] = meshgrid(  xx, xx  );
    x = reshape(X,1,P);
    y = reshape(Y,1,P);
end
pos = [x;y];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate 1D C^alpha signal
nn = 1024;          % super resolution

% generate a super resolution signal
t = 0:1/(nn-1):1;
[Y,X] = meshgrid(t,t);
% uniform C^alpha
f = gen_signal_2d(nn, alpha);      % horizontal variation

v = interp2(Y,X,f, x,y);
v = v(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform transform : split isotropic vs. axis aligned
clear options;
options.part_type = '2axis';
w1 = perform_alpert_transform_2d(v,pos,k,1,options);
err1 = l2error(w1);

clear options;
options.part_type = 'kmeans';
w2 = perform_alpert_transform_2d(v,pos,k,1,options);
err2 = l2error(w2);

% reglin
fiton = 1:floor(0.2*P);
A = (1:P)';
B = err1;
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
plot(log2(1:P), log2(err1), 'b', log2(1:P), log2(err2), 'r.');
legend('axis aligned split', 'k-means split');
axis tight;
axis([0,log2(P), -10, mmax(log2(err1))]);
title(thetitle);

% plot(log2(1:P), log2(errmin), 'k:');
plot(log2(1:P), reglin, 'k:');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform transform : 2D sliced
disp('--> Performing transform : 2D sliced');
err = [];
clear lg;
for i = 1:5
    s = 2^(i-1);
    disp(sprintf('    %d slices.', s));
    lg{i} = sprintf('%d slices', s);
    w = perform_alpert_transform_2d_sliced(v,pos,k,s);
    err = [err, l2error(w)];
end

errmin = min(err')';

% err = [err, errmin];
% lg{i+1} = 'min';

% err = [err, l2error(v)];


% reglin
fiton = 1:floor(0.2*P);
A = (1:P)';
B = errmin;
aa = log2(A(fiton));
bb = log2(B(fiton));
[coef,S] = polyfit(aa,bb,1);
disp( sprintf('--> sliced 2D Alpert transform, decreasing=%.4f', coef(1)) );
reglin = log2(A)*coef(1)+coef(2);

% title
thetitle = sprintf('\\alpha=%d, k=%d, decreasing=%.4f', alpha, k, coef(1) );

% plot
clf;
hold on;
plot(log2(1:P), log2(err));
legend(lg);
axis tight;
axis([0,log2(P), -10, max(max(log2(err(:,1:2))))]);
title(thetitle);

% plot(log2(1:P), log2(errmin), 'k:');
plot(log2(1:P), reglin, 'r--');
hold off;

hold off;