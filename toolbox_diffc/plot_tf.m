function plot_tf(T,M)

% plot_tf - plot a tensorial field.
%
%   plot_tf(T,M);
%
%   Display min/max eigenvalue and eigensystem.
%   Should be used for debug purpose, not fine-tuned display.
%
%   See also: plot_tensor_field.
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<2
    M = [];
end

[e1,e2,l1,l2] = perform_tensor_decomp(T);

clf;

subplot(2,2,1);
imagesc(l1);
title('Min eigenvalue.');
axis square;
axis off;

subplot(2,2,2);
imagesc(l2);
title('Max eigenvalue.');
axis square;
axis off;

subplot(2,2,3);
plot_vf(e1, M, 0);
title('Max eigenvector.');

subplot(2,2,4);
plot_vf(e2, M, 0);
title('Max eigenvector.');