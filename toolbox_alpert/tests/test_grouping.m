% test for recursive dichotomic grouping
%
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D sampling
n = 310;
X = rand(1,n);

part = {1:n};
clf;
for i=1:6
    clear options;
%     options.ptmax = 20;
    options.depthmax = i-1;
    options.part = part;
    options.type = '1axis';
    [part,B,G] = dichotomic_grouping(X,options);
    subplot(2,3,i);
    plot_dichotomic_partition([X;X], part);
    axis off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D sampling
n = 210;
X = rand(2,n);

part = {1:n};
clf;
for i=1:6
    clear options;
%     options.ptmax = 20;
    options.depthmax = i-1;
    options.part = part;
    options.type = '2axis';
    [part,B,G] = dichotomic_grouping(X,options);
    subplot(2,3,i);
    plot_dichotomic_partition(X, part);
    axis off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian clusters
d = 2;      % dimension
N = 300;    % points in each groups
sigma = { [1,1], [0.1,1.5], [1.5 0.1], [0.2 1.3] };
center = { [0,0] [3 3] [3 -3] [6 0]};

X = [];
for i=1:length(sigma)
    s = sigma{i};
    c = center{i};
    Xk = randn(d,N);
    for j=1:d
        Xk(j,:) = Xk(j,:)*s(j)+c(j);
    end
    X = [X, Xk];
end

n = size(X,2);
part = {1:n};
clf;
for i=1:6
    clear options;
%     options.ptmax = 20;
    options.depthmax = i-1;
    options.nb_iter = 15;
    options.part = part;
    options.type = 'kmeans';
    [part,B] = dichotomic_grouping(X,options);
    subplot(2,3,i);
    plot_dichotomic_partition(X, part);
    axis off;
end

% saveas(gcf, 'dichotomic_grouping', 'eps');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uniform sampling
n = 500;
X = rand(2,n);

clf;
part = {1:n};
for i=1:6
    clear options;
%     options.ptmax = 20;
    options.depthmax = i-1;
    options.nb_iter = 20;
    options.part = part;
    options.type = 'kmeans';
    [part,B] = dichotomic_grouping(X,options);
    subplot(2,3,i);
    plot_dichotomic_partition(X, part);
    axis off;
end

% saveas(gcf, 'dichotomic_grouping_uniform', 'eps');
