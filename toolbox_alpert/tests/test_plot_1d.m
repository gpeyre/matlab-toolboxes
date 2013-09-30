% test_plot_1d - plot 1D Alpert basis.
%
%   Copyright (c) 2004 Gabriel Peyré


k = 4;  % vanishing moments
n = 512;
x = rand(n,1);
x = sort(x);
% uncomment to have regular sampling
x = 0:1/(n-1):1;
V = build_alpert_matrix_1d(x,k);
% plot the vectors
j_list = 7;
clf;
kj = 0;
for j=j_list
kj = kj + 1;
kn = 0;
% n_list = floor( rand(1,2)*(2^(j-1)-1) )+1;
n_list = 1:4;
for i=n_list
    kn = kn+1;
    v = V(:, end-n/2^(j-1)+i );
    subplot(length(j_list), length(n_list), kn + (kj-1)*length(n_list));
    sel = 1:n/2^(4-j);
    sel = find(v~=0);
    plot(x(sel), v(sel));  
    axis tight;   
    axis square;
    if i==1
        title(sprintf('j=%d', -j));
    end
end
end

% saveas(gcf, 'alpert_basis', 'png');
disp('Press any key.');
pause;

clf;
imagesc(-abs(V));
axis square; 
set(gca, 'XTick', []);
set(gca, 'YTick', []);
colormap gray(256);
title('1D Alpert basis vectors.');
disp('Press any key.');
pause;