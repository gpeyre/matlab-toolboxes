% test mex function for haar and lifting transforms

n = 32;
Jmin = 0;
x = rand(n,1);
% x = ones(n,1);

y = perform_79_transform(x, Jmin, 1);
xx = perform_79_transform(y, Jmin, -1);
norme(x-xx)


type = '7_9';
options.use_mex = 1;
y1 = perform_lifting_transform_byname(x, Jmin, 1, type, options);
options.use_mex = 0;
y2 = perform_lifting_transform_byname(x, Jmin, 1, type, options);
norme(y1-y2)