% test_speed_2d - test the speed of the transform algorithms.
%
%   Copyright (c) 2004 Gabriel Peyré

sampling_type = 'regular';
sampling_type = 'irregular';
% position on a rectangular grid

n = 1000;
if strcmp(sampling_type, 'regular')
%    J = log2(n/k^2);
%    n = k^2*2^J;
    
    nx = sqrt(n);
    ny = sqrt(n);
    [Y,X] = meshgrid( 0:1/(nx-1):1, 0:1/(ny-1):1 );
    pos = [ reshape(X,1,n); reshape(Y,1,n) ];
    x = pos(1,:);
    y = pos(2,:);
else 
%    J = log2(n/k^2);
%    n = k^2*2^J;
    
    pos = rand(2,n);
%    [x,I] = sort(pos(1,:));
%    pos(1,:) = pos(1,I);
%    pos(2,:) = pos(2,I);
    x = pos(1,:);
    y = pos(2,:);
end

s = 1;
x = 2*(x-0.5); y = 2*(y-0.5);
v = cos(x.^2 + y.^2)';
k = 4;

% test for speed
options.part_type = '2axis';
tic;
V = build_alpert_matrix_2d(pos,k, options);
w1 = (V')*v;
disp(['2d matrix takes ' num2str(toc)]);
tic;

w2 = perform_alpert_transform_2d(v,pos,k, 1, options);
disp(['2d mex takes ' num2str(toc)]);
tic;

v1 = perform_alpert_transform_2d(w2,pos,k, -1, options);
disp(['difference between original and reconstruction (should be 0) : ' num2str(norme(v-v1))]);